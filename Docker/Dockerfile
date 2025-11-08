# Use the official R Shiny image
FROM rocker/shiny:latest

# Install any required R packages
RUN R -e "install.packages(c('tidyverse', 'shiny', 'shinydashboard'))"

# Copy your Shiny app into the image
COPY ./app /srv/shiny-server/app

# Expose port 3838
EXPOSE 3838

# Run Shiny Server
CMD ["/usr/bin/shiny-server"]

