# Use an official Ubuntu base image
FROM ubuntu:latest

# Set environment variables to prevent prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update the package list and install essential packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
        build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /app

# Copy the requirements.txt file into the container
COPY requirements.txt ./

# Install the Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy the rest of the application code into the container
COPY . ./

# Expose the default port for Streamlit
EXPOSE 8501

# Set the default command to run your Streamlit application
CMD ["streamlit", "run", "app.py"]