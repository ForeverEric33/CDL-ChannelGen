# Channel Modeling and CIR Generation Tool
This project is a **Channel Impulse Response (CIR) and Delay Channel Parameter Generator** based on the **Clustered Delay Line (CDL) model**. It is designed to simulate and generate realistic wireless channel characteristics for advanced communication systems, particularly in the context of **5G, MIMO (Multiple Input Multiple Output), and millimeter-wave (mmWave) technologies**.

The CDL model is a widely adopted channel modeling approach that captures the multi-path propagation effects, including **cluster-based delay spreads, angular spreads, and power distributions**, which are critical for accurately representing real-world wireless environments. This tool generates **time-domain CIRs** and extracts key channel parameters such as **delay spread, power delay profile (PDP), and multi-path component characteristics**.


## Key Features

- **Multi-Scenario Support**: Supports a wide range of scenarios, including:
  - Urban Macro (UMa)
  - Urban Micro (UMi)
  - Rural Macro (RMa)
  - Indoor Factory (InF)
  - Indoor Hotspot (InH)
  - And more (e.g., InF-SL, InF-DL, InF-SH, InF-DH).

- **Multi-Antenna Configurations**: Allows flexible configuration of antenna parameters for both BS and UE, including:
  - Antenna orientation (rotation angles)
  - Polarization (type and angle)
  - Directivity gain
  - Panel arrangements.

- **Dynamic Mobility Modeling**: Supports dynamic UE mobility with customizable speed and position settings.

- **CDL Model Implementation**: Faithfully implements the 3GPP-standardized CDL model, capturing multi-path propagation effects such as:
  - Cluster-based delay spreads
  - Angular spreads
  - Power distributions.
## Applications

- **Wireless System Design**: Provides a realistic channel model for testing and optimizing communication algorithms, such as channel estimation, equalization, and beamforming.
- **MIMO and mmWave Research**: Supports the study of spatial channel characteristics and multi-antenna techniques.
- **5G Simulation**: Enables the evaluation of 5G systems under various propagation conditions.
- **Hardware Testing**: Generates CIRs for validating the performance of hardware implementations (e.g., FPGAs, DSPs).

## How to Use
run main.cpp and change your own configurations.
