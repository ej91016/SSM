required_packages <- c("ggplot2", "dplyr", "tidyr","cowplot")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}


plot_ssmr <- function(mydetail = NULL, myresult, distance = "CID", mytitle = "", myrange = c(0, 0),
                      parsimony = FALSE, xlabel="", show_legend = FALSE, coverage = 0.5) {
  
  if (is.null(mydetail)) {
    mydetail <- myresult[0, , drop = FALSE]
  }
  myresult <- myresult |> dplyr::mutate(dplyr::across(-filename, ~ as.numeric(.)))
  mydetail <- mydetail |> dplyr::mutate(dplyr::across(-filename, ~ as.numeric(.)))

  lowerbound <- (1-coverage)/2
  upperbound <- 1-lowerbound
  

  if (myrange[[2]] == 0) {
    if (!is.null(mydetail)){
      myrange <- minmax(mydetail, distance, parsimony, coverage = coverage)
    }else{
      myrange <- minmax(myresult, distance, parsimony)
    }

  }
  
  type_order <- c("MP", "SSAMP", "SSA", "SSM")
  legend_order <- c("SSM", "SSA", "MP", "SSAMP")
  

  df_mean <- tibble::tibble(i = myresult$i) |>
    dplyr::mutate(
      SSM = myresult[[paste0(distance, '_SSM')]],
      SSA = myresult[[paste0(distance, '_SSA')]],
      MP = if (parsimony) myresult[[paste0(distance, '_MP')]] else NA,
      SSAMP = if (parsimony) myresult[[paste0(distance, '_SSAMP')]] else NA
    ) |>
    tidyr::pivot_longer(cols = c("SSM", "SSA", "MP", "SSAMP"), names_to = "Type", values_to = "Value") |>
    dplyr::filter(!is.na(.data$Value)) |>
    dplyr::mutate(Type = factor(Type, levels = type_order))
  

  df_raw <- tibble::tibble(i = mydetail$i) |>
    dplyr::mutate(
      SSM = mydetail[[paste0(distance, '_SSM')]],
      SSA = mydetail[[paste0(distance, '_SSA')]],
      MP = if (parsimony) mydetail[[paste0(distance, '_MP')]] else NA,
      SSAMP = if (parsimony) mydetail[[paste0(distance, '_SSAMP')]] else NA
    ) |>
    tidyr::pivot_longer(cols = c("SSM", "SSA", "MP", "SSAMP"), names_to = "Type", values_to = "Value") |>
    dplyr::filter(!is.na(.data$Value)) |>
    dplyr::mutate(Type = factor(Type, levels = type_order))
  
  shape_map <- c(MP = 17, SSAMP = 15, SSA = 16, SSM = 4)
  color_map <- c(MP = "#FFD700", SSAMP = "#007A5B", SSA = "#005F8C", SSM = "#FF1493")
  
  present_types <- intersect(legend_order, unique(df_mean$Type))
  df_mean$Type <- factor(df_mean$Type, levels = present_types)
  df_raw$Type <- factor(df_raw$Type, levels = present_types)
  
  ggplot2::ggplot() +
    ggplot2::stat_summary(data = df_raw, 
                          ggplot2::aes(x = i, y = Value, fill = Type, group = Type),
                          fun.min = function(x) quantile(x, lowerbound, na.rm = TRUE), 
                          fun.max = function(x) quantile(x, upperbound, na.rm = TRUE), 
                          geom = "ribbon", alpha = 0.2, color = NA) +
    

    ggplot2::geom_point(data = df_mean, 
                        ggplot2::aes(x = i, y = Value, color = Type, shape = Type),
                        size = 2, stroke = 0.7) +
    

    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                                labels = scales::label_number(accuracy = 0.01)) +
    ggplot2::coord_cartesian(ylim = myrange) + 
    ggplot2::labs(title = mytitle, y = distance, x = xlabel) +
    ggplot2::scale_shape_manual(values = shape_map[present_types]) +
    ggplot2::scale_color_manual(values = color_map[present_types]) +
    ggplot2::scale_fill_manual(values = color_map[present_types]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      legend.position = if (show_legend == FALSE) "none" else show_legend,
      legend.title = ggplot2::element_blank()
    )
}

plot_ssmd <- function(mydetail, myresult, distance = "CID", mytitle = "", myrange = c(0, 0),
                      parsimony = FALSE, xlabel="", show_legend = FALSE) {
  myresult <- myresult |>
    dplyr::mutate(dplyr::across(-filename, ~ as.numeric(.)))
  mydetail <- mydetail |>
    dplyr::mutate(dplyr::across(-filename, ~ as.numeric(.)))
  if (myrange[[2]] == 0) {
    myrange <- range(c(
      mydetail[[paste0(distance, "_SSA")]],
      mydetail[[paste0(distance, "_SSM")]],
      myresult[[paste0(distance, "SSA")]],
      myresult[[paste0(distance, "_SSM")]],
      if (parsimony) c(
        mydetail[[paste0(distance, "_MP")]],
        mydetail[[paste0(distance, "_SSAMP")]],
        myresult[[paste0(distance, "_MP")]],
        myresult[[paste0(distance, "_SSAMP")]]
      ) else NULL
    ), na.rm = TRUE)
  }
  df <- dplyr::bind_rows(
    tibble::tibble(i = mydetail$i, Value = mydetail[[paste0(distance, "_SSA")]], Type = "SSA_bg"),
    tibble::tibble(i = mydetail$i, Value = mydetail[[paste0(distance, "_SSM")]], Type = "SSM_bg"),
  )
  if (parsimony) {
    df <- dplyr::bind_rows(
      df,
      tibble::tibble(i = mydetail$i, Value = mydetail[[paste0(distance, "_SSAMP")]], Type = "SSAMP_bg"),
      tibble::tibble(i = mydetail$i, Value = mydetail[[paste0(distance, "_MP")]], Type = "MP_bg"),
      tibble::tibble(i = myresult$i, Value = myresult[[paste0(distance, "_SSAMP")]], Type = "SSAMP"),
      tibble::tibble(i = myresult$i, Value = myresult[[paste0(distance, "_MP")]], Type = "MP")
    )
  }
  df <- dplyr::bind_rows(
    df,
    tibble::tibble(i = myresult$i, Value = myresult[[paste0(distance, "_SSA")]], Type = "SSA"),
    tibble::tibble(i = myresult$i, Value = myresult[[paste0(distance, "_SSM")]], Type = "SSM")
  )

  color_map <- c(
    SSA     = "#005F8C",
    SSM     = "#FF1493",
    SSAMP   = "#007A5B",
    MP      = "#FFD700",
    
    SSA_bg   = "#005F8C",
    SSM_bg   = "#FF1493",
    SSAMP_bg = "#007A5B",
    MP_bg    = "#FFD700"
  )
  
  shape_map <- c(
    SSM = 4, SSA = 16, MP = 17, SSAMP = 15,
    SSM_bg = 4, SSA_bg = 16, MP_bg = 17, SSAMP_bg = 15
  )
  
  alpha_map <- c(
    SSM = 1, SSA = 1, MP = 1, SSAMP = 1,
    SSM_bg = 0.05, SSA_bg = 0.05, MP_bg = 0.05, SSAMP_bg = 0.05
  )
  type_order <- c("SSA_bg", "SSAMP_bg", "MP_bg", "SSM_bg", "SSA", "SSAMP", "MP", "SSM")
  legend_order <- c("SSM","SSA","MP","SSAMP","SSM_bg","SSA_bg","MP_bg","SSAMP_bg")
  present_types <- intersect(legend_order, unique(df$Type))
  df$Type <- factor(df$Type, levels = present_types)
  
  ggplot2::ggplot(df, ggplot2::aes(x = i, y = Value, color = Type, shape = Type, alpha = Type)) +
    ggplot2::geom_point(size = 2,stroke=0.7) +
    ggplot2::scale_shape_manual(values = shape_map[present_types]) +
    ggplot2::scale_color_manual(values = color_map[present_types]) +
    ggplot2::scale_alpha_manual(values = alpha_map[present_types]) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::labs(title = mytitle, y = distance, x = xlabel) +
    ggplot2::ylim(myrange) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      legend.position = if (show_legend == FALSE) "none" else show_legend,
      legend.title = ggplot2::element_blank()
    )
}

plot_legend <- function(parsimony = FALSE) {
  library(ggplot2)
  

  dummy_df <- data.frame(
    i = 1:4,
    Value = 1:4,
    Type = factor(
      c("SSM", "SSA", "MP", "SSAMP")[1:(2 + 2 * as.integer(parsimony))],
      levels = c("SSM", "SSA", "MP", "SSAMP")
    )
  )
  
  shape_map <- c(SSM = 4, SSA = 16, MP = 17, SSAMP = 15)
  color_map <- c(    SSA     = "#005F8C",
                     SSM     = "#FF1493",
                     SSAMP   = "#007A5B",
                     MP      = "#FFD700"
              )
  

  p <- ggplot(dummy_df, aes(x = i, y = Value, color = Type, shape = Type)) +
    geom_point(size = 3,stroke=0.7) +
    scale_shape_manual(values = shape_map) +
    scale_color_manual(values = color_map) +
    labs(shape = NULL, color = NULL) +
    theme_void() +
    theme(legend.position = "right")
  cowplot::get_legend(p)
}

minmax <- function(mydata, distance = 'CID', parsimony = FALSE, coverage = 1) {
  lower_q <- (1 - coverage) / 2
  upper_q <- 1 - lower_q
  
  cols <- c(paste0(distance, '_SSM'), paste0(distance, '_SSA'))
  if (parsimony) {
    cols <- c(cols, paste0(distance, '_MP'), paste0(distance, '_SSAMP'))
  }
  
  per_step_bounds <- mydata |>
    dplyr::group_by(i) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(cols), 
                                   list(low = ~quantile(as.numeric(.), lower_q, na.rm=TRUE),
                                        high = ~quantile(as.numeric(.), upper_q, na.rm=TRUE))))
  

  global_min <- min(per_step_bounds |> dplyr::select(dplyr::contains("_low")), na.rm = TRUE)
  global_max <- max(per_step_bounds |> dplyr::select(dplyr::contains("_high")), na.rm = TRUE)
  
  return(c(0.95 * global_min, 1.05 * global_max))
}
