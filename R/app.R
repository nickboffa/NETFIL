library(shiny)
library(leaflet)
library(plotly)
library(dplyr)
library(sf)

source("R/prep_output_for_shiny.R")

# ── Pre-load datasets ──────────────────────────────────────────────────────────
DATASETS <- list(
  "Many"      = prep_csv_for_shiny("Many",      "many_3mda_x15.csv"), # many_test for 5x 0MDA
  "Raster660" = prep_csv_for_shiny("Raster660", "r660_3mda_x15.csv") # raster660_test for 5x 0MDA
)

# ── User-adjustable parameters ─────────────────────────────────────────────────
PREV_MAX    <- 7
MS_PER_STEP <- 267

# ── Per-model colours ──────────────────────────────────────────────────────────
MODEL_LINE_COL   <- c("Many" = "#7ecfff", "Raster660" = "#ffb347")
MODEL_RIBBON_COL <- c("Many" = "rgba(126,207,255,0.12)", "Raster660" = "rgba(255,179,71,0.12)")

# ── Pre-compute mean + 90% range across sims for every model (done once) ───────
df_all_overall <- bind_rows(lapply(names(DATASETS), function(nm) {
  DATASETS[[nm]] %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    group_by(SimulationNumber, Time) %>%
    summarise(Mf_sim = mean(Mf_prev, na.rm = TRUE), .groups = "drop") %>%
    mutate(Time_round = round(Time, 2)) %>%
    group_by(Time_round) %>%
    summarise(
      Mf_mean = mean(Mf_sim,            na.rm = TRUE),
      Mf_lo   = quantile(Mf_sim, 0.05,  na.rm = TRUE),
      Mf_hi   = quantile(Mf_sim, 0.95,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Model = nm)
}))

Y_MAX_ALL   <- max(df_all_overall$Mf_hi, na.rm = TRUE) * 1.1
# Vline trace index (0-based): 2 traces per model (ribbon+line) + 1 current-sim line + vline
VLINE_TRACE <- as.integer(2L * length(DATASETS) + 1L)

# ── Helpers ────────────────────────────────────────────────────────────────────
group_label <- function(g) {
  if (suppressWarnings(!is.na(as.numeric(g)))) "" else paste0("<b>Group:</b> ", g, "<br>")
}

pal <- colorNumeric(
  palette  = c("#2c7bb6", "#ffffbf", "#d7191c"),
  domain   = c(0, PREV_MAX),
  na.color = "grey"
)
clamp_prev <- function(x) pmin(pmax(x, 0), PREV_MAX)

initial_time_steps <- DATASETS[[1]] %>%
  mutate(Time_round = round(Time, 2)) %>%
  pull(Time_round) %>% unique() %>% sort()

initial_sims <- sort(unique(DATASETS[[1]]$SimulationNumber))

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      * { box-sizing: border-box; }
      html, body { margin: 0; padding: 0; height: 100%; overflow: hidden;
                   font-family: sans-serif; background: #1a1a2e; color: #eee; }

      #outer-grid {
        display: grid;
        grid-template-columns: 1fr 2fr;
        grid-template-rows: 42vh 58vh;
        height: 100vh; width: 100vw; gap: 0;
      }
      #cell-stats { grid-column: 1; grid-row: 1; padding: 20px;
                    display: flex; flex-direction: column; justify-content: center; }
      #cell-plot  { grid-column: 2; grid-row: 1; padding: 10px 16px 0 0; min-width: 0; }
      #cell-map   { grid-column: 1 / span 2; grid-row: 2; position: relative; }

      #stats-box {
        background: rgba(126,207,255,0.06); border: 1px solid #7ecfff33;
        border-radius: 14px; padding: 24px 28px; height: 100%;
        display: flex; flex-direction: column; justify-content: center;
      }
      .stat-overall-label { font-size: 11px; color: #7ecfff; letter-spacing: .12em;
                             text-transform: uppercase; margin-bottom: 4px; }
      .stat-overall { font-size: 52px; font-weight: 800; color: #fff;
                      line-height: 1; margin-bottom: 4px; }
      .stat-overall-unit { font-size: 22px; color: #aac; font-weight: 400; }
      hr.stat-divider { border: none; border-top: 1px solid #7ecfff22; margin: 16px 0; }
      .stat-row { display: flex; justify-content: space-between; align-items: center;
                  font-size: 12px; color: #889; margin-top: 10px; }
      .stat-val { color: #7ecfff; font-weight: 700; font-size: 16px; }

      #cell-plot .plotly { height: 100% !important; }
      #map { height: 100%; width: 100%; }
      .leaflet-control-zoom { display: none !important; }

      #controls {
        position: absolute; bottom: 18px; left: 50%;
        transform: translateX(-50%);
        background: rgba(20,20,40,0.88); backdrop-filter: blur(8px);
        padding: 12px 22px 10px; border-radius: 12px; z-index: 1000;
        width: min(760px, 92%); box-shadow: 0 4px 24px rgba(0,0,0,0.5);
      }
      #model-row {
        display: flex; justify-content: center; align-items: center;
        gap: 14px; margin-bottom: 6px; flex-wrap: wrap;
      }
      #model-row label { color: #7ecfff; font-size: 12px; letter-spacing: .08em;
                         text-transform: uppercase; margin: 0; white-space: nowrap; }
      #model_select {
        background: rgba(126,207,255,0.1); color: #eee;
        border: 1px solid #7ecfff44; border-radius: 6px;
        padding: 3px 8px; font-size: 13px;
      }
      /* numericInput: strip Shiny wrapper, style the <input> directly */
      #sim_select { margin: 0; }
      #sim_select input[type=number] {
        background: rgba(126,207,255,0.1) !important; color: #eee !important;
        border: 1px solid #7ecfff44 !important; border-radius: 6px !important;
        padding: 3px 8px !important; font-size: 13px !important;
        height: 28px !important; width: 72px !important;
        -moz-appearance: textfield !important;
      }
      #sim_select input[type=number]::-webkit-inner-spin-button,
      #sim_select input[type=number]::-webkit-outer-spin-button { -webkit-appearance: none; }
      .form-group { margin-bottom: 0 !important; }

      #time_display { text-align: center; font-size: 16px; font-weight: bold;
                      color: #7ecfff; margin-bottom: 4px; letter-spacing: .08em; }
      .btn-play { background: #7ecfff; color: #1a1a2e; border: none;
                  padding: 6px 16px; border-radius: 6px; font-weight: bold;
                  cursor: pointer; font-size: 13px; white-space: nowrap; }
      .btn-play:hover { background: #a8e0ff; }
      .btn-reset { background: rgba(126,207,255,0.15); color: #7ecfff;
                   border: 1px solid #7ecfff44; padding: 6px 10px;
                   border-radius: 6px; font-size: 15px; cursor: pointer; line-height: 1; }
      .btn-reset:hover { background: #7ecfff; color: #1a1a2e; }
      .slider-row { display: flex; align-items: center; gap: 8px; margin-top: 4px; }
      .slider-row input[type=range] { flex: 1; accent-color: #7ecfff; }
      .irs--shiny .irs-min, .irs--shiny .irs-max,
      .irs--shiny .irs-single { display: none !important; }
      .speed-btns { display: flex; gap: 6px; justify-content: center; margin-top: 6px; }
      .btn-speed { background: rgba(126,207,255,0.15); color: #7ecfff;
                   border: 1px solid #7ecfff44; padding: 3px 12px;
                   border-radius: 5px; font-size: 12px; cursor: pointer; }
      .btn-speed:hover, .btn-speed.active { background: #7ecfff; color: #1a1a2e; }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('setSpeedActive', function(id) {
        document.querySelectorAll('.btn-speed').forEach(b => b.classList.remove('active'));
        document.getElementById(id).classList.add('active');
      });
    "))
  ),

  div(id = "outer-grid",

      div(id = "cell-stats",
          div(id = "stats-box",
              div(class = "stat-overall-label", "Overall Mf prevalence"),
              div(class = "stat-overall",
                  textOutput("stat_overall", inline = TRUE),
                  tags$span(class = "stat-overall-unit", "%")
              ),
              tags$hr(class = "stat-divider"),
              div(class = "stat-row",
                  span("Highest group"),
                  span(class = "stat-val", textOutput("stat_max", inline = TRUE))
              ),
              div(class = "stat-row",
                  span("Lowest group"),
                  span(class = "stat-val", textOutput("stat_min", inline = TRUE))
              )
          )
      ),

      div(id = "cell-plot",
          plotlyOutput("prev_plot", width = "100%", height = "100%")
      ),

      div(id = "cell-map",
          leafletOutput("map", width = "100%", height = "100%"),
          absolutePanel(
            id = "controls", bottom = 18, left = "50%",
            div(id = "model-row",
                tags$label(`for` = "model_select", "Model output"),
                selectInput("model_select", label = NULL,
                            choices = names(DATASETS),
                            selected = names(DATASETS)[1],
                            width = "160px"),
                tags$label(`for` = "sim_select",
                           textOutput("sim_label", inline = TRUE)),
                numericInput("sim_select", label = NULL,
                             value = min(initial_sims),
                             min   = min(initial_sims),
                             max   = max(initial_sims),
                             step  = 1, width = "72px")
            ),
            div(id = "time_display", textOutput("time_label", inline = TRUE)),
            div(class = "slider-row",
                actionButton("play_btn",  "▶ Play", class = "btn-play"),
                actionButton("reset_btn", "⏮",     class = "btn-reset"),
                sliderInput("time_slider", label = NULL,
                            min = 1, max = length(initial_time_steps), value = 1,
                            step = 1, animate = FALSE, width = "100%", ticks = FALSE)
            ),
            div(class = "speed-btns",
                actionButton("spd1", "1×", class = "btn-speed active"),
                actionButton("spd2", "2×", class = "btn-speed"),
                actionButton("spd5", "5×", class = "btn-speed")
            )
          )
      )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Clean base data for selected model ────────────────────────────────────
  dataset <- reactive({
    DATASETS[[input$model_select]] %>% filter(!is.na(lat), !is.na(lon))
  })

  time_steps <- reactive({
    dataset() %>%
      mutate(Time_round = round(Time, 2)) %>%
      pull(Time_round) %>% unique() %>% sort()
  })

  # ── Validated sim number ───────────────────────────────────────────────────
  sim_val <- reactive({
    v <- input$sim_select
    req(!is.null(v), !is.na(v))
    sims <- sort(unique(dataset()$SimulationNumber))
    as.integer(max(min(round(v), max(sims)), min(sims)))
  })

  output$sim_label <- renderText({
    sims <- sort(unique(dataset()$SimulationNumber))
    paste0("Sim (0–", max(sims), ")")
  })

  # ── Map data: selected simulation ─────────────────────────────────────────
  df_map <- reactive({
    dataset() %>%
      filter(SimulationNumber == sim_val()) %>%
      group_by(Group, Time, lat, lon) %>%
      summarise(mf = mean(mf, na.rm = TRUE), pop = mean(pop, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(Mf_prev = 100 * mf / pop, Time_round = round(Time, 2))
  })

  # ── Current-sim line: per-time mean across groups ──────────────────────────
  df_current_sim <- reactive({
    dataset() %>%
      filter(SimulationNumber == sim_val()) %>%
      group_by(Time) %>%
      summarise(Mf_sim = mean(Mf_prev, na.rm = TRUE), .groups = "drop") %>%
      mutate(Time_round = round(Time, 2))
  })

  # ── Reset when model changes ───────────────────────────────────────────────
  observeEvent(input$model_select, {
    sims <- sort(unique(dataset()$SimulationNumber))
    updateNumericInput(session, "sim_select",
                       value = min(sims), min = min(sims), max = max(sims))
    playing(FALSE)
    updateActionButton(session, "play_btn", label = "▶ Play")
    updateSliderInput(session, "time_slider", max = length(time_steps()), value = 1)
  })

  # ── Speed buttons ──────────────────────────────────────────────────────────
  speed_mult <- reactiveVal(1)
  observe_speed <- function(id, mult) {
    observeEvent(input[[id]], {
      speed_mult(mult)
      session$sendCustomMessage("setSpeedActive", id)
    })
  }
  observe_speed("spd1", 1); observe_speed("spd2", 2); observe_speed("spd5", 5)

  playing <- reactiveVal(FALSE)
  observeEvent(input$reset_btn, {
    playing(FALSE)
    updateActionButton(session, "play_btn", label = "▶ Play")
    updateSliderInput(session, "time_slider", value = 1)
  })

  current_data <- reactive({
    t <- time_steps()[input$time_slider]
    df_map() %>% filter(Time_round == t)
  })
  current_time <- reactive({ time_steps()[input$time_slider] })
  output$time_label <- renderText({ paste("Year:", round(current_time(), 1)) })

  # ── Stats (from selected simulation) ──────────────────────────────────────
  output$stat_overall <- renderText({
    round(mean(current_data()$Mf_prev, na.rm = TRUE), 2)
  })
  output$stat_max <- renderText({
    paste0(round(max(current_data()$Mf_prev, na.rm = TRUE), 2), "%")
  })
  output$stat_min <- renderText({
    paste0(round(min(current_data()$Mf_prev, na.rm = TRUE), 2), "%")
  })

  # ── Plotly ─────────────────────────────────────────────────────────────────
  # Trace layout (0-indexed):
  #   0,1   ribbon+line model 1
  #   2,3   ribbon+line model 2  (repeats for N models)
  #   2N    current-sim line
  #   2N+1  vline  ← VLINE_TRACE
  output$prev_plot <- renderPlotly({
    ts     <- time_steps()
    df_sim <- df_current_sim()
    init_t <- ts[isolate(input$time_slider)]

    p <- plot_ly(source = "prev_plot")

    for (nm in names(DATASETS)) {
      dfo <- filter(df_all_overall, Model == nm)
      p <- p %>%
        add_ribbons(data = dfo,
                    x = ~Time_round, ymin = ~Mf_lo, ymax = ~Mf_hi,
                    fillcolor = MODEL_RIBBON_COL[nm],
                    line = list(color = "transparent"),
                    hoverinfo = "none", showlegend = FALSE) %>%
        add_lines(data = dfo,
                  x = ~Time_round, y = ~Mf_mean,
                  name = nm,
                  line = list(color = MODEL_LINE_COL[nm], width = 2),
                  hovertemplate = paste0(nm, " – %{x:.1f}: %{y:.2f}%<extra></extra>"))
    }

    p %>%
      add_lines(data = df_sim,
                x = ~Time_round, y = ~Mf_sim,
                name = paste0(input$model_select, " sim ", sim_val()),
                line = list(color = "rgba(255,255,255,0.6)", width = 1.5, dash = "dash"),
                hovertemplate = "%{x:.1f}: %{y:.2f}%<extra></extra>") %>%
      add_segments(x = init_t, xend = init_t, y = 0, yend = Y_MAX_ALL,
                   line = list(color = "#ff7e7e", width = 2, dash = "dot"),
                   inherit = FALSE, showlegend = FALSE) %>%
      layout(
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor  = "rgba(0,0,0,0)",
        font   = list(color = "#aac", family = "sans-serif", size = 11),
        xaxis  = list(title = "Year", gridcolor = "#7ecfff11",
                      tickfont = list(color = "#aac"), zeroline = FALSE),
        yaxis  = list(title = "Mf prevalence (%)", gridcolor = "#7ecfff11",
                      tickfont = list(color = "#aac"), zeroline = FALSE,
                      rangemode = "tozero"),
        legend = list(font = list(color = "#aac", size = 10),
                      bgcolor = "rgba(20,20,40,0.6)",
                      bordercolor = "#7ecfff22", borderwidth = 1,
                      x = 0.01, y = 0.99, xanchor = "left", yanchor = "top"),
        annotations = list(list(
          text = "Shaded area: 90% range across simulations",
          x = 1, y = 0, xref = "paper", yref = "paper",
          xanchor = "right", yanchor = "bottom",
          showarrow = FALSE,
          font = list(size = 9, color = "#55667788")
        )),
        margin = list(t = 10, r = 10, b = 44, l = 50),
        showlegend = TRUE
      ) %>%
      config(displayModeBar = FALSE)
  })

  observeEvent(event_data("plotly_click", source = "prev_plot"), {
    click <- event_data("plotly_click", source = "prev_plot")
    if (!is.null(click)) {
      nearest_idx <- which.min(abs(time_steps() - click$x[1]))
      updateSliderInput(session, "time_slider", value = nearest_idx)
    }
  })

  observeEvent(input$time_slider, {
    t <- current_time()
    plotlyProxy("prev_plot", session) %>%
      plotlyProxyInvoke("restyle",
                        list(x = list(c(t, t)), y = list(c(0, Y_MAX_ALL))),
                        VLINE_TRACE)
  })

  # ── Leaflet base map ───────────────────────────────────────────────────────
  output$map <- renderLeaflet({
    dm <- df_map()
    leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
      addProviderTiles(providers$CartoDB.DarkMatter) %>%
      fitBounds(
        lng1 = min(dm$lon, na.rm = TRUE), lat1 = min(dm$lat, na.rm = TRUE) - 0.1,
        lng2 = max(dm$lon, na.rm = TRUE), lat2 = max(dm$lat, na.rm = TRUE) + 0.01
      ) %>%
      addLegend(position = "topright", pal = pal,
                values = seq(0, PREV_MAX, length.out = 100),
                title  = "Mf prevalence (%)", opacity = 0.85)
  })

  # ── Markers ────────────────────────────────────────────────────────────────
  observe({
    cd <- current_data()
    sqrt_pop <- sqrt(pmax(cd$pop, 1))
    radii    <- 2 * sqrt_pop / mean(sqrt_pop, na.rm = TRUE)  # average radius = 2

    labels <- lapply(seq_len(nrow(cd)), function(i) {
      htmltools::HTML(paste0(
        group_label(cd$Group[i]),
        "<b>Mf+:</b> ", round(cd$Mf_prev[i], 2), "% (",
        round(cd$mf[i]), "/", round(cd$pop[i]), ")"
      ))
    })

    leafletProxy("map") %>%
      addCircleMarkers(
        data = cd, lng = ~lon, lat = ~lat, layerId = ~Group,
        radius      = radii,
        color       = ~pal(clamp_prev(Mf_prev)),
        fillColor   = ~pal(clamp_prev(Mf_prev)),
        fillOpacity = 0.85, stroke = TRUE, weight = 1, opacity = 0.6,
        group = "points", label = labels,
        labelOptions = labelOptions(
          style = list("background" = "rgba(20,20,40,0.85)", "color" = "#eee",
                       "border" = "none", "border-radius" = "6px",
                       "padding" = "6px 10px", "font-size" = "13px"),
          direction = "auto", textOnly = FALSE
        )
      )
  })

  # ── Playback ───────────────────────────────────────────────────────────────
  play_timer <- reactiveTimer(MS_PER_STEP)

  observeEvent(input$play_btn, {
    playing(!playing())
    updateActionButton(session, "play_btn",
                       label = if (playing()) "⏸ Pause" else "▶ Play")
  })

  observe({
    play_timer()
    if (!playing()) return()
    ts        <- isolate(time_steps())
    current   <- isolate(input$time_slider)
    mult      <- isolate(speed_mult())
    next_step <- min(current + mult, length(ts))
    updateSliderInput(session, "time_slider", value = next_step)
    if (next_step >= length(ts)) {
      playing(FALSE)
      updateActionButton(session, "play_btn", label = "▶ Play")
    }
  })
}

shinyApp(ui, server)
