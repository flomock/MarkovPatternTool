# Use an official Python runtime as a base image
#FROM python:2-slim
FROM python:2

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
ADD . /app

#RUN apk --update add --virtual build-dependencies gcc


#RUN apt-get update && apt-get install -y apt-transport-https

#RUN apt-get update \
#    && apt-get install -y --no-install-recommends apt-transport-https gcc and-build-dependencies \
#    && rm -rf /var/lib/apt/lists/* \
#    && pip install cryptography \
#    && apt-get purge -y --auto-remove gcc and-build-dependencies

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME World

# Run app.py when the container launches
CMD ["python", "FractalMatrix.py"]
