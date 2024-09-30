# Use the rocker/tidyverse base image
FROM rocker/tidyverse:latest

# Install system dependencies for your package (if any)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev

# Install necessary R packages
RUN R -e "install.packages(c('devtools', 'remotes'))"

# Copy your R package source code into the container
COPY . /admixtools

# Set the working directory
WORKDIR /admixtools

# Install your R package
RUN R -e "devtools::install_local('.')"

# Verify the installation
RUN R -e "if (!requireNamespace('admixtools', quietly = TRUE)) { stop('admixtools package installation failed') }"

CMD ["R"]

