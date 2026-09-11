using PlotlyJS

craig_lab_template = PlotlyJS.Template(
    layout = attr(
        autotypenumbers =  "strict",
        colorway = [
            "#8b5058",
            "#525b7a",
            "#ec9f79",
            "#a781a2",
            "#9fc2cc",
            "#d17466"
        ],
        font = attr(
            color =  "#2a3f5f",
            size =  18,
            family =  "Tahoma"
        ),
        hovermode = "closest",
        hoverlabel = attr(
            align = "left"
        ),
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        polar= attr(
            bgcolor= "white",
            angularaxis= attr(
                gridcolor = "#EBF0F8",
                linecolor = "#EBF0F8",
                ticks = ""
            ),
            radialaxis = attr(
                gridcolor = "#EBF0F8",
                linecolor = "#EBF0F8",
                ticks = ""
            )
        ),
        ternary = attr(
            bgcolor = "white",
            aaxis = attr(
                gridcolor = "#DFE8F3",
                linecolor = "#A2B1C6",
                ticks = ""
            ),
            baxis = attr(
                gridcolor = "#DFE8F3",
                linecolor = "#A2B1C6",
                ticks = ""
            ),
            caxis = attr(
                gridcolor = "#DFE8F3",
                linecolor = "#A2B1C6",
                ticks = ""
            )
        ),
        coloraxis = attr(
            colorbar = attr(
                outlinewidth = 0,
                ticks = ""
            )
        ),
        colorscale = attr(
            sequential = [
                [0, "#05172C"],
                [0.07142857142857142, "#0E2C49"],
                [0.14285714285714285, "#224163"],
                [0.21428571428571427, "#455779"],
                [0.2857142857142857, "#5B5C78"],
                [0.35714285714285715, "#715F74"],
                [0.42857142857142855, "#88616E"],
                [0.5, "#A06268"],
                [0.5714285714285714, "#BB6461"],
                [0.6428571428571429, "#D57262"],
                [0.7142857142857143, "#E58165"],
                [0.7857142857142857, "#E79A75"],
                [0.8571428571428571, "#E8B087"],
                [0.9285714285714286, "#EDD1AA"],
                [1, "#FCF2D5"]
            ],
            sequentialminus = [
                [0, "#FCF2D5"],
                [0.07142857142857142, "#EDD1AA"],
                [0.14285714285714285, "#E8B087"],
                [0.21428571428571427, "#E79A75"],
                [0.2857142857142857, "#E58165"],
                [0.35714285714285715, "#D57262"],
                [0.42857142857142855, "#BB6461"],
                [0.5, "#A06268"],
                [0.5714285714285714, "#88616E"],
                [0.6428571428571429, "#715F74"],
                [0.7142857142857143, "#5B5C78"],
                [0.7857142857142857, "#455779"],
                [0.8571428571428571, "#224163"],
                [0.9285714285714286, "#0E2C49"],
                [1, "#05172C"]
            ],
            diverging = [
                [0.0, "#292f40"],
                [0.07142857142857142, "#3e475d"],
                [0.14285714285714285, "#556177"],
                [0.21428571428571427, "#6d7d91"],
                [0.2857142857142857, "#8699ab"],
                [0.35714285714285715, "#a1b7c4"],
                [0.42857142857142855, "#c0d5dd"],
                [0.5, "#f6ebe3"],
                [0.5714285714285714, "#eec6bd"],
                [0.6428571428571429, "#dda49a"],
                [0.7142857142857143, "#c6857a"],
                [0.7857142857142857, "#ab675e"],
                [0.8571428571428571, "#8f4c45"],
                [0.9285714285714286, "#72332e"],
                [1.0, "#551b19"]
            ]
        ),
        xaxis = attr(
            gridcolor = "#EBF0F8",
            linecolor = "#C9C3C1",
            tickcolor = "#C9C3C1",
            ticks = "inside",
            title = attr(
                standoff = 15,
                font = attr(size = 18)
            ),
            zerolinecolor = "#EBF0F8",
            automargin = true,
            zerolinewidth = 2,
            showline = true,
            mirror = true,
            showgrid = false
        ),
        yaxis = attr(
            gridcolor = "#EBF0F8",
            linecolor = "#C9C3C1",
            tickcolor = "#C9C3C1",
            ticks = "inside",
            title = attr(
                standoff = 15,
                font = attr(size = 18)
            ),
            zerolinecolor = "#EBF0F8",
            automargin = true,
            zerolinewidth = 2,
            showline = true,
            mirror = true,
            showgrid = false
        ),
        scene = attr(
            xaxis = attr(
                backgroundcolor = "white",
                gridcolor = "#DFE8F3",
                linecolor = "#EBF0F8",
                showbackground = true,
                ticks = "",
                zerolinecolor = "#EBF0F8",
                gridwidth = 2
            ),
            yaxis = attr(
                backgroundcolor = "white",
                gridcolor = "#DFE8F3",
                linecolor = "#EBF0F8",
                showbackground = true,
                ticks = "",
                zerolinecolor = "#EBF0F8",
                gridwidth = 2
            ),
            zaxis = attr(
                backgroundcolor = "white",
                gridcolor = "#DFE8F3",
                linecolor = "#EBF0F8",
                showbackground = true,
                ticks = "",
                zerolinecolor = "#EBF0F8",
                gridwidth = 2
            )
        ),
        shapedefaults = attr(
            line = attr(
                color = "#525b7a"
            )
        ),
        annotationdefaults = attr(
            arrowcolor = "#525b7a",
            arrowhead = 0,
            arrowwidth = 1,
            font = attr(size = 24)
        ),
        geo = attr(
            bgcolor = "white",
            landcolor = "white",
            subunitcolor = "#C8D4E3",
            showland = true,
            showlakes = true,
            lakecolor = "white"
        ),
        title = attr(
            x = 0.05,
            font = attr(size = 20)
        ),
        mapbox = attr(
            style = "light"
        )
    ),
    data = Dict([
        (:histogram2dcontour,
            [ attr(
                type = "histogram2dcontour",
                colorbar =attr(
                    outlinewidth = 0,
                    ticks = ""),
                colorscale = [
                    [0, "#05172C"],
                    [0.07142857142857142, "#0E2C49"],
                    [0.14285714285714285, "#224163"],
                    [0.21428571428571427, "#455779"],
                    [0.2857142857142857, "#5B5C78"],
                    [0.35714285714285715, "#715F74"],
                    [0.42857142857142855, "#88616E"],
                    [0.5, "#A06268"],
                    [0.5714285714285714, "#BB6461"],
                    [0.6428571428571429, "#D57262"],
                    [0.7142857142857143, "#E58165"],
                    [0.7857142857142857, "#E79A75"],
                    [0.8571428571428571, "#E8B087"],
                    [0.9285714285714286, "#EDD1AA"],
                    [1, "#FCF2D5"]
                ]
            )
            ]
        ),
        (:choropleth,
            [ attr(
                type = "choropleth",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = "")
                )
            ]
        ),
        (:histogram2d,
            [ attr(
                type = "histogram2d",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = ""),
                colorscale = [
                    [0, "#05172C"],
                    [0.07142857142857142, "#0E2C49"],
                    [0.14285714285714285, "#224163"],
                    [0.21428571428571427, "#455779"],
                    [0.2857142857142857, "#5B5C78"],
                    [0.35714285714285715, "#715F74"],
                    [0.42857142857142855, "#88616E"],
                    [0.5, "#A06268"],
                    [0.5714285714285714, "#BB6461"],
                    [0.6428571428571429, "#D57262"],
                    [0.7142857142857143, "#E58165"],
                    [0.7857142857142857, "#E79A75"],
                    [0.8571428571428571, "#E8B087"],
                    [0.9285714285714286, "#EDD1AA"],
                    [1, "#FCF2D5"]
                ])
            ]
        ),
        (:heatmap,
            [   attr(
                type = "heatmap",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = "")
                ,
                colorscale = [
                    [0, "#05172C"],
                    [0.07142857142857142, "#0E2C49"],
                    [0.14285714285714285, "#224163"],
                    [0.21428571428571427, "#455779"],
                    [0.2857142857142857, "#5B5C78"],
                    [0.35714285714285715, "#715F74"],
                    [0.42857142857142855, "#88616E"],
                    [0.5, "#A06268"],
                    [0.5714285714285714, "#BB6461"],
                    [0.6428571428571429, "#D57262"],
                    [0.7142857142857143, "#E58165"],
                    [0.7857142857142857, "#E79A75"],
                    [0.8571428571428571, "#E8B087"],
                    [0.9285714285714286, "#EDD1AA"],
                    [1, "#FCF2D5"]
                ]
            )
            ]
        ),
        (:heatmapgl, [
            attr(
                type = "heatmapgl",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = ""
                ),
                colorscale = [
                    [0, "#05172C"],
                    [0.07142857142857142, "#0E2C49"],
                    [0.14285714285714285, "#224163"],
                    [0.21428571428571427, "#455779"],
                    [0.2857142857142857, "#5B5C78"],
                    [0.35714285714285715, "#715F74"],
                    [0.42857142857142855, "#88616E"],
                    [0.5, "#A06268"],
                    [0.5714285714285714, "#BB6461"],
                    [0.6428571428571429, "#D57262"],
                    [0.7142857142857143, "#E58165"],
                    [0.7857142857142857, "#E79A75"],
                    [0.8571428571428571, "#E8B087"],
                    [0.9285714285714286, "#EDD1AA"],
                    [1, "#FCF2D5"]
                ]
            )
        ]),
        (:contourcarpet, [
            attr(
                type = "contourcarpet",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = ""
                )
            )
        ]),
        (:contour, [
            attr(
                type = "contour",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = ""
                ),
                colorscale = [
                    [0, "#05172C"],
                    [0.07142857142857142, "#0E2C49"],
                    [0.14285714285714285, "#224163"],
                    [0.21428571428571427, "#455779"],
                    [0.2857142857142857, "#5B5C78"],
                    [0.35714285714285715, "#715F74"],
                    [0.42857142857142855, "#88616E"],
                    [0.5, "#A06268"],
                    [0.5714285714285714, "#BB6461"],
                    [0.6428571428571429, "#D57262"],
                    [0.7142857142857143, "#E58165"],
                    [0.7857142857142857, "#E79A75"],
                    [0.8571428571428571, "#E8B087"],
                    [0.9285714285714286, "#EDD1AA"],
                    [1, "#FCF2D5"]
                ]
            )
        ]),
        (:surface, [
            attr(
                type = "surface",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = ""
                ),
                colorscale = [
                    [0, "#05172C"],
                    [0.07142857142857142, "#0E2C49"],
                    [0.14285714285714285, "#224163"],
                    [0.21428571428571427, "#455779"],
                    [0.2857142857142857, "#5B5C78"],
                    [0.35714285714285715, "#715F74"],
                    [0.42857142857142855, "#88616E"],
                    [0.5, "#A06268"],
                    [0.5714285714285714, "#BB6461"],
                    [0.6428571428571429, "#D57262"],
                    [0.7142857142857143, "#E58165"],
                    [0.7857142857142857, "#E79A75"],
                    [0.8571428571428571, "#E8B087"],
                    [0.9285714285714286, "#EDD1AA"],
                    [1, "#FCF2D5"]
                ]
            )
        ]),
        (:mesh3d, [
            attr(
                type = "mesh3d",
                colorbar = attr(
                    outlinewidth = 0,
                    ticks = ""
                )
            )
        ]),
        (:scatter,
            [ attr(
                type = "scatter",
                marker = attr(
                    size = 8,
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                ),
                line = attr(width = 3)
                )
            ]
        ),
        (:parcoords, [
            attr(
                type = "parcoords",
                line = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:scatterpolargl, [
            attr(
                type = "scatterpolargl",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:bar, [
            attr(
                error_x = attr(
                    color = "#2a3f5f"
                ),
                error_y = attr(
                    color = "#2a3f5f"
                ),
                marker = attr(
                    line = attr(
                        color = "white",
                        width = 0.5
                    )
                ),
                type = "bar"
            )
        ]),
        (:scattergeo, [
            attr(
                type = "scattergeo",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:scatterpolar, [
            attr(
                type = "scatterpolar",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:histogram, [
            attr(
                type = "histogram",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:scattergl, [
            attr(
                type = "scattergl",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:scatter3d, [
            attr(
                type = "scatter3d",
                line = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                ),
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:scattermapbox, [
            attr(
                type = "scattermapbox",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:scatterternary, [
            attr(
                type = "scatterternary",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:scattercarpet, [
            attr(
                type = "scattercarpet",
                marker = attr(
                    colorbar = attr(
                        outlinewidth = 0,
                        ticks = ""
                    )
                )
            )
        ]),
        (:carpet, [
            attr(
                aaxis = attr(
                    endlinecolor = "#2a3f5f",
                    gridcolor = "#C8D4E3",
                    linecolor = "#C8D4E3",
                    minorgridcolor = "#C8D4E3",
                    startlinecolor = "#2a3f5f"
                ),
                baxis = attr(
                    endlinecolor = "#2a3f5f",
                    gridcolor = "#C8D4E3",
                    linecolor = "#C8D4E3",
                    minorgridcolor = "#C8D4E3",
                    startlinecolor = "#2a3f5f"
                ),
                type = "carpet"
            )
        ]),
        (:table, [
            attr(
                cells = attr(
                    fill = attr(
                        color = "#EBF0F8"
                    ),
                    line = attr(
                        color = "white"
                    )
                ),
                header = attr(
                    fill = attr(
                        color = "#C8D4E3"
                    ),
                    line = attr(
                        color = "white"
                    )
                ),
                type = "table"
            )
        ]),
        (:barpolar, [
            attr(
                marker = attr(
                    line = attr(
                        color = "white",
                        width = 0.5
                    )
                ),
                type = "barpolar"
            )
        ]),
        (:pie, [
            attr(
                automargin = true,
                type = "pie"
            )
        ])
    ])
)
