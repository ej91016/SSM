required_packages <- c("ggplot2", "dplyr", "tidyr","cowplot")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}

plot_ssmr <- function(myresult, distance = "CID", mytitle = "", myrange = c(0, 0),
                      parsimony = FALSE, xlabel="", show_legend = FALSE) {
  myresult <- myresult |>
    dplyr::mutate(dplyr::across(-filename, ~ as.numeric(.)))

  if (myrange[[2]] == 0) {
    myrange <- minmax(myresult, distance, parsimony)
  }
  type_order <- c("MP", "SSAMP", "SSA", "SSM")
  legend_order <- c("SSM", "SSA", "MP", "SSAMP")
  df <- tibble::tibble(i = myresult$i) |>
    dplyr::mutate(
      SSM = myresult[[distance]],
      SSA = myresult[[paste0(distance, '_p')]],
      MP = if (parsimony) myresult[[paste0(distance, '_pp')]] else NA,
      SSAMP = if (parsimony) myresult[[paste0(distance, '_pps')]] else NA
    ) |>
    tidyr::pivot_longer(cols = c("SSM", "SSA", "MP", "SSAMP"), names_to = "Type", values_to = "Value") |>
    dplyr::filter(!is.na(.data$Value)) |>
    dplyr::mutate(Type = factor(Type, levels = type_order)) |>
    dplyr::arrange(Type)
  

  shape_map <- c(MP = 17, SSAMP = 15, SSA = 16, SSM = 4)
  color_map <- c(MP = "#FFD700", SSAMP = "#007A5B", SSA = "#005F8C", SSM = "#FF1493")

  present_types <- intersect(legend_order, unique(df$Type))
  df$Type <- factor(df$Type, levels = present_types)
  
  
  ggplot2::ggplot(df, ggplot2::aes(x = i, y = Value, color = Type, shape = Type)) +
    ggplot2::geom_point(size = 2,stroke=0.7) +
    ggplot2::scale_shape_manual(values = shape_map[present_types]) +
    ggplot2::scale_color_manual(values = color_map[present_types]) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::scale_y_continuous(
      limits = myrange,
      breaks = scales::pretty_breaks(n = 5),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    ggplot2::labs(title = mytitle, y = distance, x = xlabel) +
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
      mydetail[[paste0(distance, "_p")]],
      mydetail[[distance]],
      myresult[[paste0(distance, "_p")]],
      myresult[[distance]],
      if (parsimony) c(
        mydetail[[paste0(distance, "_pp")]],
        mydetail[[paste0(distance, "_pps")]],
        myresult[[paste0(distance, "_pp")]],
        myresult[[paste0(distance, "_pps")]]
      ) else NULL
    ), na.rm = TRUE)
  }
  df <- dplyr::bind_rows(
    tibble::tibble(i = mydetail$i, Value = mydetail[[paste0(distance, "_p")]], Type = "SSA_bg"),
    tibble::tibble(i = mydetail$i, Value = mydetail[[distance]], Type = "SSM_bg"),
  )
  if (parsimony) {
    df <- dplyr::bind_rows(
      df,
      tibble::tibble(i = mydetail$i, Value = mydetail[[paste0(distance, "_pps")]], Type = "SSAMP_bg"),
      tibble::tibble(i = mydetail$i, Value = mydetail[[paste0(distance, "_pp")]], Type = "MP_bg"),
      tibble::tibble(i = myresult$i, Value = myresult[[paste0(distance, "_pps")]], Type = "SSAMP"),
      tibble::tibble(i = myresult$i, Value = myresult[[paste0(distance, "_pp")]], Type = "MP")
    )
  }
  df <- dplyr::bind_rows(
    df,
    tibble::tibble(i = myresult$i, Value = myresult[[paste0(distance, "_p")]], Type = "SSA"),
    tibble::tibble(i = myresult$i, Value = myresult[[distance]], Type = "SSM")
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
  
  # Dummy data
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
  
  # Build dummy plot to extract the legend
  p <- ggplot(dummy_df, aes(x = i, y = Value, color = Type, shape = Type)) +
    geom_point(size = 3,stroke=0.7) +
    scale_shape_manual(values = shape_map) +
    scale_color_manual(values = color_map) +
    labs(shape = NULL, color = NULL) +
    theme_void() +
    theme(legend.position = "right")
  cowplot::get_legend(p)
}

minmax <- function(mydata, distance = 'CID', parsimony = FALSE) {
  if (parsimony) {
    return(c(
      0.95 * min(
        as.numeric(mydata[[distance]]),
        as.numeric(mydata[[paste0(distance, '_p')]]),
        as.numeric(mydata[[paste0(distance, '_pp')]]),
        as.numeric(mydata[[paste0(distance, '_pps')]])
      ),
      1.05 * max(
        as.numeric(mydata[[distance]]),
        as.numeric(mydata[[paste0(distance, '_p')]]),
        as.numeric(mydata[[paste0(distance, '_pp')]]),
        as.numeric(mydata[[paste0(distance, '_pps')]])
      )
    ))
  } else {
    return(c(
      0.95 * min(
        as.numeric(mydata[[distance]]),
        as.numeric(mydata[[paste0(distance, '_p')]])
      ),
      1.05 * max(
        as.numeric(mydata[[distance]]),
        as.numeric(mydata[[paste0(distance, '_p')]])
      )
    ))
  }
}
