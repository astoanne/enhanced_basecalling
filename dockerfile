# Use an official lightweight base image
FROM alpine:latest

# Set the working directory
WORKDIR /app

# Create a script to print "Hello, World!"
RUN echo 'echo "Hello, World!"' > hello.sh && \
    chmod +x hello.sh

# Set the default command to execute the script
CMD ["./hello.sh"]