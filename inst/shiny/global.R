library(sass)

sass(
    sass_file("styles/main.scss"),
    options = sass_options(output_style = "compressed"),
    cache = NULL,
    output = "www/css/dali.min.css"
)
