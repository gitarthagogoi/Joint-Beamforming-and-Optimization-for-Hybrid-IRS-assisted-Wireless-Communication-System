# Hybrid IRS MU-MISO Simulation (Work-in-Progress)

This repository contains an early-stage MATLAB simulation framework for a **Hybrid Intelligent Reflecting Surface (Active + Passive IRS)** assisted MU-MISO wireless communication system.  
The goal is to build a modular testbed for studying hybrid IRS behavior, system power consumption, and performance trade-offs in next-generation wireless networks.

## 🎯 Project Objective
To model and analyze:
- Hybrid IRS operation (active + passive elements)
- Spectral Efficiency (SE) and Energy Efficiency (EE)
- Power consumption including BS + IRS hardware
- Interaction between BS precoding and IRS element control

The current focus is establishing the core simulation pipeline and verifying consistent behavior across channel realizations.

---

## ✅ Current Features
- Channel modeling (BS–IRS, IRS–User, BS–User)
- Hybrid IRS parameter structure (gain + phase)
- Total power consumption model
- SE and EE computation
- Multi-user MU-MISO support
- Initial framework for hybrid precoder + IRS optimization

---

## 🚧 Work Remaining
### **Documentation**
- Model assumptions  
- System architecture overview  
- Function-level documentation  
- Usage instructions  

### **Hybrid Precoding Implementation**
Upcoming major development:
- Hybrid BS precoding module  
- Joint optimization of **W** (BS precoder) and **Θ** (Hybrid IRS matrix)  
- Verification experiments to validate the hybrid scheme  

This will complete the full hybrid IRS evaluation phase.

---

## 📌 Status
**Under active development.**  
Core functions work, but documentation and hybrid precoder verification are still pending.

---
