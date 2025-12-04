Hybrid IRS MU-MISO Simulation (Work-in-Progress)

This repository contains an early-stage MATLAB simulation framework for a Hybrid Intelligent Reflecting Surface (Active + Passive IRS) assisted MU-MISO wireless communication system.
The project aims to explore hybrid IRS modeling, channel characterization, and power-aware performance evaluation for next-generation wireless networks.

Goal of the Project

The primary objective is to build a flexible testbed that can evaluate:

Hybrid IRS behavior (active + passive elements)

SE and EE under hardware-aware power consumption

Interaction between BS precoding and IRS element control

Influence of amplifier constraints in active IRS panels

This repo currently focuses on establishing the core simulation pipeline, ensuring that channel modeling, power accounting, and metric evaluation behave consistently across multiple realizations.

Current Features

Channel modeling for BS–IRS, IRS–User, and BS–User links

Hybrid IRS parameter setup (active amplification + passive phase shifts)

Total power consumption modeling (transmit + IRS hardware)

Spectral Efficiency (SE) and Energy Efficiency (EE) evaluation

Support for multi-user MU-MISO configurations

Early structure for joint precoder/IRS optimization

Work Remaining

This is an active development project. The following components are planned next:

🔧 Documentation

Clear explanation of model assumptions

Hybrid IRS architecture overview

Function-level descriptions

Usage guide and simulation flow diagram

🧠 Hybrid Precoding Implementation (Upcoming)

The next major task is implementing and verifying:

Hybrid BS precoding scheme

Joint optimization of W (BS precoder) and Θ (IRS hybrid phase–gain matrix)

Verification experiments to validate the hybrid design

This phase will complete the full hybrid IRS evaluation pipeline.

Status

🚧 Under Development — Code is functional but not fully verified.
Expect updates as optimization modules and documentation are added.
