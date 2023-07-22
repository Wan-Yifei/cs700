# Use Ubuntu 20.04 LTS as the base image
FROM ubuntu:20.04

# Update package lists and install necessary packages
RUN apt-get update && apt-get install -y \
    python3.8 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Set the default Python version
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1

# Set the working directory
WORKDIR /app

# Copy application code into the container
COPY . /app
COPY ./bin /bin

# Install dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Set the entrypoint file
#ENTRYPOINT ["/app/entrypoint.sh"]
