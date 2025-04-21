required_packages <- c("ape", "phytools", "TreeDist", "phangorn", 
                       "processx", "foreach", "doParallel", "Matrix", "readr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}

clustergen <- function(npar,rpar){
  num_cores = min(npar,parallel::detectCores()%/%rpar - 1)
  print('code will be excuted in parallel')
  print(paste0(num_cores,'*',rpar,' threads in use'))
  cl = parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  return(cl)
}

char_gen <-function(tree_num,state_num,set_num,mytrees,rate,ntip,char_num,
                    falsespace=NULL,char_total,fullDir,parsimony=FALSE){
  for (t in 1:tree_num){
    mytree=mytrees[[t]]
    myname=paste0(mytree$tip.label," ")
    mymodels=list()
    for (k in 1:length(state_num)){
      mymodels[[k]]=mkQmatrix(n=state_num[[k]],r=rate)
    }
    for (i in 1:set_num){
      mymatrix=data.frame(matrix(nrow=ntip,ncol=0))
      for (k in 1:length(state_num)){
        if (char_num[k]!=0){
          myparts=phytools::sim.Mk(mytree,mymodels[[k]],nsim=char_num[[k]])
          if(k==1){
            mymatrix=data.frame(myparts)
          }else{
            mymatrix=cbind(mymatrix,myparts)
          }
        }
      }
      if (!is.null(falsespace)){
        dum=replicate(ntip, 0:(falsespace-1))
        mymatrix=cbind(mymatrix,t(dum))
      }
      fname=file.path(fullDir,paste0('physim_',t,'_',i,'.phy'))
      title=paste(' ',ntip,' ',char_total,"\n",sep="")
      cat(title,file=fname)
      write.table(mymatrix,file=fname,row.names = myname,col.names = FALSE
                  ,quote = FALSE,sep="",append=T)
      if (!is.null(falsespace)){
        stnum_mod=falsespace
      }else{
        stnum_mod=max(state_num)
      }
      if (parsimony){
        makepaup(fullDir=fullDir,pt=t,pi=i,pnt=ntip,pctotal=char_total,pst=stnum_mod,
                 pmymatrix=mymatrix, pname=myname)
      }
    }
  }
}

par_char_gen <- function(tree_num,state_num,set_num,mytrees,rate,ntip,char_num,
                         falsespace, char_total,fullDir, parsimony,cl) {
  parallel::clusterExport(cl, c("mkQmatrix", "makepaup"))
  `%dopar%` <- foreach::`%dopar%`
  tryCatch({
    foreach::foreach(t = 1:tree_num, .packages = c("ape", "phytools","Matrix")) %dopar% {
      mytree <- mytrees[[t]]
      myname <- paste0(mytree$tip.label, " ")
      mymodels <- list()
      for (k in 1:length(state_num)) {
        mymodels[[k]] <- mkQmatrix(n = state_num[[k]], r = rate)
      }
      for (i in 1:set_num) {
        mymatrix <- data.frame(matrix(nrow = ntip, ncol = 0))
        for (k in 1:length(state_num)) {
          if (char_num[k] != 0) {
            myparts <- phytools::sim.Mk(mytree, mymodels[[k]], nsim = char_num[[k]])
            if (k == 1) {
              mymatrix <- data.frame(myparts)
            } else {
              mymatrix <- cbind(mymatrix, myparts)
            }
          }
        }
        if (!is.null(falsespace)) {
          dum <- replicate(ntip, 0:(falsespace - 1))
          mymatrix <- cbind(mymatrix, t(dum))
        }
        # Write matrix to file
        fname <- file.path(fullDir, paste0('physim_', t, '_', i, '.phy'))
        cat(fname)
        title <- paste(' ', ntip, ' ', char_total, "\n", sep = "")
        cat(title, file = fname)
        write.table(mymatrix,file=fname,row.names=myname,col.names=FALSE,
                    quote=FALSE,sep="",append=TRUE)
        
        if (!is.null(falsespace)) {
          stnum_mod <- falsespace
        } else {
          stnum_mod <- max(state_num)
        }
        if (parsimony) {
          makepaup(fullDir=fullDir,pt=t,pi=i,pnt=ntip,pctotal=char_total,
                   pst=stnum_mod,pmymatrix=mymatrix,pname=myname)
        }
      }
    }    
  },error=function(e){
    message("Error in function par_char_gen: ",e)
    })
}

makepaup <- function(fullDir,pt,pi,pnt,pctotal,pst,pmymatrix,pname){
  f2name=paste0('paupsim_',pt,'_',pi,'.nex')
  paupsize=paste0('dimensions ntax=',pnt,' nchar=',pctotal,';')
  paupchar=paste0('FORMAT DATATYPE = STANDARD GAP = - MISSING = ? SYMBOLS = " ',
                  paste0(c(0:(pst-1)),collapse=' '),'";')
  pauplog=paste0('log file = paup_',pt,'_',pi,'.txt replace;')
  paupsave=paste0('contree / format=newick treeFile=', fullDir,'/paup_physim_',
                  pt,'_',pi,'.tre replace;')
  cat(c('#nexus','begin data;',paupsize,paupchar,'matrix\n'),
      file=file.path(fullDir,f2name),sep='\n')
  write.table(pmymatrix,file = file.path(fullDir,f2name), row.names = pname, col.names = FALSE
              ,quote = FALSE,sep="",append=T)
  cat(c('\n;','end;\n','begin paup;\n\n',pauplog,
        'set warntree=no warnreset=no warnTSave=no increase=auto autoInc=100;',
        'hsearch reconLimit=20 rearrLimit=10000000 allSwap=yes limitPerRep=yes;',
        'describeTrees 1;',paupsave,'quit;'), file=file.path(fullDir,f2name),sep='\n',append=T)
}

#on average 1 substitution per site
mkQmatrix <- function(n=2,rate=1){
  qmatrix=matrix(c(rep(1/(n-1),n^2)),n,dimnames=list(0:(n-1),0:(n-1)))
  diag(qmatrix)<- (-1)
  qmatrix=qmatrix*rate
  return(qmatrix)
}

raxml_run <- function(threads, model, mfile, myDir) {
  for (filename in mfile) {
    treename <- paste0(tools::file_path_sans_ext(filename), ".tre")
    # Basic RAxML run
    processx::run(
      "raxmlHPC-PTHREADS-SSE3",
      args = c(
        "-m", "MULTICAT",
        "-T", threads,
        "-K", "MK",
        "-p", "12345",
        "-s", file.path(myDir, filename),
        "-V",
        "-n", treename
      ),
      wd = myDir,
      stdout = file.path(myDir,"raxmlrun_check.txt"),
      stderr = file.path(myDir,"raxmlrun_error.txt")
    )
    
    # Partitioned run
    processx::run(
      "raxmlHPC-PTHREADS-SSE3",
      args = c(
        "-m", "MULTICAT",
        "-T", threads,
        "-K", "MK",
        "-p", "12345",
        "-s", file.path(myDir, filename),
        "-V",
        "-q", model,
        "-n", paste0("part_", treename)
      ),
      wd = myDir,
      stdout = file.path(myDir,"raxmlrun_check.txt"),
      stderr = file.path(myDir,"raxmlrun_error.txt")
    )
  }
}

par_raxml_run <- function(threads, model,mfile,myDir,cl) {
  # Parallel loop using foreach to process each file
  `%dopar%` <- foreach::`%dopar%`
  tryCatch({
    foreach::foreach(filename = mfile, .packages = c("processx")) %dopar% {
      log_file = paste0("log_", basename(filename), ".txt")
      treename <- paste0(tools::file_path_sans_ext(filename), ".tre")
      # Basic RAxML run
      processx::run(
        "raxmlHPC-PTHREADS-SSE3",
        args = c(
          "-m", "MULTICAT",
          "-T", threads,
          "-K", "MK",
          "-p", "12345",
          "-s", file.path(myDir, filename),
          "-V",
          "-n", treename
        ),
        wd = myDir,
        stdout = file.path(myDir,"raxmlrun_check.txt"),
        stderr = file.path(myDir,"raxmlrun_error.txt")
      )
      
      # Partitioned run
      processx::run(
        "raxmlHPC-PTHREADS-SSE3",
        args = c(
          "-m", "MULTICAT",
          "-T", threads,
          "-K", "MK",
          "-p", "12345",
          "-s", file.path(myDir, filename),
          "-V",
          "-q", model,
          "-n", paste0("part_", treename)
        ),
        wd = myDir,
        stdout = file.path(myDir,"raxmlrun_check.txt"),
        stderr = file.path(myDir,"raxmlrun_error.txt")
      )
    }
  },error=function(e){
    message("Error in function par_raxml_run: ",e)
  })
}

paup_run <- function(mfile,myDir){
  for (filename in mfile){
    processx::run('cmd', 
                  args = c('/c', 'paup',filename),
                  wd = myDir,
                  stdout = file.path(myDir,'pauprun_check.txt'),
                  stderr = file.path(myDir,'pauprun_error.txt'))
  }
}

par_paup_run <- function(mfile,myDir,cl) {
  # Parallel loop using foreach to process each file
  `%dopar%` <- foreach::`%dopar%`
  tryCatch({
    foreach::foreach(filename = mfile, .packages = c("processx")) %dopar% {
      setwd(myDir)
      processx::run('cmd', 
                    args = c('/c', 'paup',filename),
                    wd = myDir,
                    stdout = 'pauprun_check.txt',
                    stderr = 'pauprun_error.txt')
    }
  },error=function(e){
    message("Error in function par_paup_run: ",e)
  })
}

tree_sum <- function(fullDir, tree_num,set_num,parsimony,rate=1){
  sim_tree=ape::read.tree(file.path(fullDir,'sim.trees'))
  for (i in 1:length(sim_tree)){
    sim_tree[[i]]$edge.length=sim_tree[[i]]$edge.length*rate
  }
  for (j in 1:tree_num){
    comp_name=file.path(fullDir,paste0('comp_physim_',j,'.trees'))
    ape::write.tree(sim_tree[j],file=comp_name)
    part_comp_name=file.path(fullDir,paste0('part_comp_physim_',j,'.trees'))
    ape::write.tree(sim_tree[j],file=part_comp_name)
    if (parsimony){
      paup_comp_name=file.path(fullDir,paste0('paup_comp_physim_',j,'.trees'))
      ape::write.tree(sim_tree[j],file=paup_comp_name)
    }
    for (i in 1:set_num){
      raxml_name=paste0('RAxML_bestTree.physim_',j,'_',i,'.tre')
      temp_string=readr::read_file(file.path(fullDir,raxml_name))
      cat(temp_string,file=comp_name,append=TRUE)
      part_raxml_name=paste0('RAxML_bestTree.part_physim_',j,'_',i,'.tre')
      part_temp_string=readr::read_file(file.path(fullDir,part_raxml_name))
      cat(part_temp_string,file=part_comp_name,append=TRUE)
      if (parsimony){
        paup_name=paste0('paup_physim_',j,'_',i,'.tre')
        paup_temp_string=readr::read_file(file.path(fullDir,paup_name))
        cat(paup_temp_string,file=paup_comp_name,append=TRUE)
      }
    }
  }
}

comp_sum <- function(fullDir,tree_num,set_num,mysetting,parsimony){
  single=c()
  single_PID=c()
  single_CID=c()
  part=c()
  part_PID=c()
  part_CID=c()
  if(parsimony){
    paup=c()
    paup_PID=c()
    paup_CID=c()
  }
  for (j in 1:tree_num){
    sumtree=paste0('comp_physim_',j,'.trees')
    comptrees=ape::read.tree(file.path(fullDir,sumtree))
    single_RF=TreeDist::RobinsonFoulds(comptrees)
    single_PIDRF=TreeDist::InfoRobinsonFoulds(comptrees,
                                    normalize=TreeDist::SplitwiseInfo(comptrees[[1]]))
    single_CIDRF=TreeDist::ClusteringInfoDist(comptrees,
                                    normalize=TreeDist::ClusteringEntropy(comptrees[[1]]))
    single=append(single,mean(single_RF[1:set_num]))
    single_PID=append(single_PID,mean(single_PIDRF[1:set_num]))
    single_CID=append(single_CID,mean(single_CIDRF[1:set_num]))
    
    part_sumtree=paste0('part_comp_physim_',j,'.trees')
    part_comptrees=ape::read.tree(file.path(fullDir,part_sumtree))
    part_RF=TreeDist::RobinsonFoulds(part_comptrees)
    part_PIDRF=TreeDist::InfoRobinsonFoulds(part_comptrees,
                                  normalize=TreeDist::SplitwiseInfo(comptrees[[1]]))
    part_CIDRF=TreeDist::ClusteringInfoDist(part_comptrees,
                                  normalize=TreeDist::ClusteringEntropy(comptrees[[1]]))
    part=append(part,mean(part_RF[1:set_num]))
    part_PID=append(part_PID,mean(part_PIDRF[1:set_num]))
    part_CID=append(part_CID,mean(part_CIDRF[1:set_num]))
    
    if(parsimony){
      paup_sumtree=paste0('paup_comp_physim_',j,'.trees')
      paup_comptrees=ape::read.tree(file.path(fullDir,paup_sumtree))
      paup_RF=TreeDist::RobinsonFoulds(paup_comptrees)
      paup_PIDRF=TreeDist::InfoRobinsonFoulds(paup_comptrees,
                                    normalize=TreeDist::SplitwiseInfo(comptrees[[1]]))
      paup_CIDRF=TreeDist::ClusteringInfoDist(paup_comptrees,
                                    normalize=TreeDist::ClusteringEntropy(comptrees[[1]]))
      paup=append(paup,mean(paup_RF[1:set_num]))
      paup_PID=append(paup_PID,mean(paup_PIDRF[1:set_num]))
      paup_CID=append(paup_CID,mean(paup_CIDRF[1:set_num]))
      dists=data.frame('RF'=single,'PID'=single_PID,'CID'=single_CID,
                       'RF_p'=part,'PID_p'=part_PID,'CID_p'=part_CID,
                       'RF_pp'=paup,'PID_pp'=paup_PID,'CID_pp'=paup_CID)
    }else{  
      dists=data.frame('RF'=single,'PID'=single_PID,'CID'=single_CID,
                       'RF_p'=part,'PID_p'=part_PID,'CID_p'=part_CID)
    }
  }
  write.csv(dists,file.path(fullDir,paste0(mysetting,'_result.csv')),row.names = FALSE)
  return(dists)
}

getresult <- function(parsimony=FALSE){
  if(parsimony){
    return(data.frame('RFD'=numeric(),'PID'=numeric(),'CID'=numeric(),
                      'RFD_p'=numeric(),'PID_p'=numeric(),'CID_p'=numeric(),
                      'RFD_pp'=numeric(),'PID_pp'=numeric(),'CID_pp'=numeric(),
                      'i'=numeric(),'filename'=character()))
  }else{
    return(data.frame('RFD'=numeric(),'PID'=numeric(),'CID'=numeric(),
                      'RFD_p'=numeric(),'PID_p'=numeric(),'CID_p'=numeric(),
                      'i'=numeric(),'filename'=character()))
  }
}

compare_file <-function(file_1,file_2,options='CID',adjust=1){
  file1=ape::read.tree(file_1)
  file2=ape::read.tree(file_2)
  reftree=file1[[1]]
  reftree$edge.length=reftree$edge.length*adjust
  mylist=c()
  mylist2=c()
  mylist3=c()
  for (i in 2:length(file1)){
    if (options=='CID'){
      dist=TreeDist::ClusteringInfoDist(file1[[i]],reftree,
                                        normalize=TreeDist::ClusteringEntropy(reftree))
      mylist=append(mylist,dist)
      dist2=TreeDist::ClusteringInfoDist(file2[[i]],reftree,
                                         normalize=TreeDist::ClusteringEntropy(reftree))
      mylist2=append(mylist2,dist2)
      dist3=TreeDist::ClusteringInfoDist(file1[[i]],file2[[i]],
                                         normalize=TreeDist::ClusteringEntropy(reftree))
      mylist3=append(mylist3,dist3)
    }
    if (options=='BSD'){
      dist=phangorn::KF.dist(file1[[i]],reftree)
      mylist=append(mylist,dist)
      dist2=phangorn::KF.dist(file2[[i]],reftree)
      mylist2=append(mylist2,dist2)
      dist3=phangorn::KF.dist(file1[[i]],file2[[i]])
      mylist3=append(mylist3,dist3)
    }
    if (options=='wRFD'){
      dist=phangorn::wRF.dist(file1[[i]],reftree)
      mylist=append(mylist,dist)
      dist2=phangorn::wRF.dist(file2[[i]],reftree)
      mylist2=append(mylist2,dist2)
      dist3=phangorn::wRF.dist(file1[[i]],file2[[i]])
      mylist3=append(mylist3,dist3)
    }
  }
  return(list(mylist,mylist2,mylist3))
}

distcal_summary<-function(dm,index,filename,result,detail){
  result[index,]=c(round(colMeans(dm),5),index,basename(filename))
  detail[(1+(index-1)*nrow(dm)):(index*nrow(dm)),1:(ncol(detail)-2)]=round(dm,5)
  detail$i[(1+(index-1)*nrow(dm)):(index*nrow(dm))]=index
  detail$filename[(1+(index-1)*nrow(dm)):(index*nrow(dm))]=basename(filename)
  return(list(result=result,detail=detail))
}

build_tree<-function(ntip,tree_num,birth,death){
  trees=c()
  while (length(trees)<tree_num){
    tt=phytools::pbtree(b=birth,d=death,n=ntip,scale=1)
    if (length(tt$tip.label)==ntip){
      if (length(trees)==0){
        trees=phytools::as.multiPhylo(tt)
      }else{
        trees=c(trees,tt)
      }
    }
  }
  return(trees)
}

datasetup<-function(falsespace,state_num,char_num,fullDir){
  if(!is.null(falsespace)){
    state_num_mod=c(state_num,falsespace)
    char_num_mod=c(char_num,falsespace)
  }else{
    state_num_mod=state_num
    char_num_mod=char_num
  }
  model_name='sim'
  for (i in 1:length(state_num)){
    model_name=paste(model_name,state_num[i],char_num[i],sep='_')
  }
  if(!is.null(falsespace)){
    model_name=paste0(model_name,'_f',falsespace)
  }
  model_name=paste0(model_name,'.models')
  
  for (i in 1:length(state_num_mod)){
    if(i==1 && char_num_mod[i]!=0){
      line=paste0('MULTI, m',i,'=','1-',char_num_mod[1],'\n')
      cat(line,file=file.path(fullDir,model_name))
    }else{
      if(char_num_mod[i]!=0){
        line=paste0('MULTI, m',i,'=',sum(char_num_mod[1:i-1])+1,'-',
                    sum(char_num_mod[1:i]),'\n')
        cat(line,file=file.path(fullDir,model_name),append=TRUE)
      }
    }
  }
  return(model_name)
}
setting<-function(ntip,birth,death,state_num,char_num,falsespace,rate){
  mysetting=paste0('trial_',ntip,'_',round(death/birth,2),'d')
  
  for (i in 1:length(state_num)){
    mysetting=paste(mysetting,state_num[i],char_num[i],sep='_')
  }
  char_total=sum(char_num)
  if(!is.null(falsespace)){
    char_total=char_total+falsespace
    mysetting=paste0(mysetting,'_f',falsespace)
  }
  mysetting=paste0(mysetting,'_r',rate)
  return(list(char_total,mysetting))
}

matrix_gen <- function(par,cl, ...) {
  if (par > 0) {
    par_char_gen(cl=cl,...)
  } else {
    char_gen(...)
  }
}

phylo_analysis <- function(par, rthreads, mymodel, fullDir,parsimony,cl) {
  raxml_file = list.files(path=fullDir,pattern = 'phy$')
  if (par>1) {
    par_raxml_run(threads=rthreads,model=mymodel,mfile=raxml_file,myDir=fullDir,cl=cl)
  } else {
    raxml_run(threads=rthreads,model=mymodel,mfile=raxml_file,myDir=fullDir)
  }
  if (parsimony) {
    paup_file = list.files(path=fullDir,pattern = 'nex$')
    if (par>0) {
      par_paup_run(paup_file,fullDir,cl=cl)
    } else {
      paup_run(paup_file,fullDir)
    }
  }
}

run_summary <- function(mainDir = getwd(), parsimony = FALSE) {
  resultframe <- getresult(parsimony=parsimony)
  detailframe <- getresult(parsimony=parsimony)
  if (Sys.which("fd") != "") {
    # Use fd if available (faster)
    rawlist <- processx::run(
      "fd",
      args = c("trial_.*result\\.csv", mainDir, "--type=file"),
      echo = FALSE,
      error_on_status = FALSE
    )$stdout
    distfile_list <- unlist(strsplit(trimws(rawlist), "\n"))
  } else {

    warning("'fd' not found - using list.files(). Install fd for better performance: https://github.com/sharkdp/fd")
    distfile_list <- list.files(
      path = mainDir,
      pattern = glob2rx("trial_*result.csv"),
      full.names = TRUE,
      recursive = TRUE
    )
  }
  

  if (length(distfile_list) == 0) {
    stop("No matching files found in: ", mainDir)
  }
  
  # Sort by creation time
  distfile_info <- file.info(distfile_list)
  if (any(is.na(distfile_info$ctime))) {
    warning("Some files lack creation times - using alphabetical order instead.")
    distfile_list <- sort(distfile_list)
  } else {
    distfile_list <- distfile_list[order(distfile_info$ctime)]
  }
  
  cat("Summarizing", length(distfile_list), "files:\n")
  cat(paste0(basename(distfile_list), collapse = "\n"), "\n")
  
  # Process files
  for (i in seq_along(distfile_list)) {
    dist_data <- read.csv(distfile_list[i])
    update <- distcal_summary(
      dm = dist_data,
      index = i,
      filename = distfile_list[i],
      result = resultframe,
      detail = detailframe
    )
    resultframe <- update$result
    detailframe <- update$detail
  }
  
  return(list(result = resultframe, detail = detailframe))
}

## ---------------------------------------------------------------------------------------------------------------------------------------------------
#this function generate one analysis for a set of the parameters
distcal <- function(state_num, char_num, rate, falsespace = NULL, ntip, birth, death,
                    tree_num, set_num, mainDir = getwd(), parsimony = FALSE, 
                    rthreads = 1, par = 0) {
  cl <- NULL
  tryCatch({

    # Generate trees and setup
    if (par > 0) {cl <- clustergen(npar = par, rpar = rthreads)
    } else {
      cat(rthreads, "threads in use")
    }
    
    mytrees <- build_tree(ntip = ntip, tree_num = tree_num, birth = birth, death = death)
    myspec <- setting(ntip, birth, death, state_num, char_num, falsespace, rate)
    
    output_dir <- file.path(mainDir, myspec[[2]])
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    } else {
      warning("Output directory already exists: ", output_dir)
    }
    
    sim_tree_path <- file.path(output_dir, "sim.trees")
    ape::write.tree(mytrees, file = sim_tree_path)
    
    # Model setup and analysis
    mymodel <- datasetup(falsespace, state_num, char_num, output_dir)
    
    matrix_gen(tree_num = tree_num, state_num = state_num, set_num = set_num,
               mytrees = mytrees, rate = rate, ntip = ntip, char_num = char_num,
               falsespace = falsespace, char_total = myspec[[1]], 
               fullDir = output_dir, parsimony = parsimony, par = par, cl = cl)
    
    phylo_analysis(par = par, rthreads = rthreads, mymodel = mymodel,
                   fullDir = output_dir, parsimony = parsimony, cl = cl)
    
    # Process results
    tree_sum(fullDir = output_dir, tree_num = tree_num, set_num = set_num,
             parsimony = parsimony, rate = rate)
    
    result <- comp_sum(fullDir = output_dir, tree_num = tree_num, set_num = set_num,
                       mysetting = myspec[[2]], parsimony = parsimony)
    
    return(list(distance = result, 
                fname = file.path(output_dir, paste0(myspec[[2]], "_result.csv"))))
    
  }, error = function(e) {
    message("Error occurred in distcal: ", e)
    return(NULL)
  }, finally = {
    if (!is.null(cl)) parallel::stopCluster(cl)
  })
}
