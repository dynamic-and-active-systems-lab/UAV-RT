## Download Instructions

Both the UAV-RT and the Pixstream folder contain programs designed for the UAV-RT project. These programs are designed to be downloaded on-board the UAV's companion computer. Pixtream contains a custom python script and xml file necessary to create a custom GNU Radio block. This block records incoming radio data, matches it to UAV data, and stores the information locally. The UAV-RT folder contains a series of custom scripts, designed to set up a UDP connection and two-way communication between a custom MATLAB GUI and the UAV. These scripts start and stop GNU Radio based on incoming messages. While the scripts can be used on their own, please refernece the [software installation instructions](https://uavrt.nau.edu/index.php/docs/radiotelem/companion-computer/) for a complete setup. 

Note:
  - must have flight controller (Pixhawk) plugged into companion computer on startup
