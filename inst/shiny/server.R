function(input, output) {

    vals <- reactiveValues(data = .GlobalEnv$.data.object.VDJ)

    dataUploadModal <- function(failed = F) {
        modalDialog(
            fileInput("file", "Choose an rds file to load", accept = ".rds"),
            if (failed) {
                div("Invalid file!")
            },

            footer = tagList(
                actionButton("load", "Load")
            ),
            easyClose = F
        )
    }

    if (is.null(isolate(vals$data)) || !isValidSeuratObject(isolate(vals$data))) {
        showModal(dataUploadModal())
    }

    observeEvent(input$load, {
        data <- readRDS(input$file$datapath)

        if (isValidSeuratObject(data)) {
            vals$data <- readRDS(input$file$datapath)
            removeModal()
        } else {
            showModal(dataUploadModal(failed = T))
        }
    })
}
