# 🧱 Finite Element Code — Dynamic Behaviour Analysis

This project provides a **finite element code** to compute the **dynamic behaviour of beam structures**, including:

- ✅ **Modal response** (natural frequencies and mode shapes)  
- ✅ **Forced response** to **harmonic excitations** (future parts)

The implementation is based on beam elements and currently supports basic lumped masses and clamped boundary conditions.  

---

## 📌 Part 1 — Modal Response Computation

### ⚙️ 1. Geometry Definition (`geometry.py`)

Define the structure's geometry by editing the following variables:

- **`nL`** — List of nodes  
  ```python
  nL[i] = [X_coordinate, Y_coordinate, Z_coordinate]
  ```
- **`eL`** — List of the base elements (linkage of your nodes) index numbering is 1-based here
  ```python
  eL[i] = [start node index of element i, end node index of element i] 
  ```
- **`etLL`** — List of type of the elements (related to their section properties ) 
  ```python
   etL[i] = type of element i
  ```
---

### 🧮 2. System Parameters (`constants.py`)

Define material, section, and discretization parameters:

- **`elemPerBar`** — Number of finite elements per beam  
  → Same number for all beams (no local refinement yet).

- **`lumpedMass`** — Total lumped mass of the structure.  

- **`desiredFreqNb`** — Number of **natural frequencies** to compute  
  → By default, the **first frequencies** are computed.

- **Material & Section Properties**  
  → Currently supports **one material** and **two section types**. This can be easily extended.

---
#### ▶️ To run the modal analysis:
```bash
python ./MECA0029_Group16_1.py
```
---

### 🧱 3. Meshing (`mesh.py`)

Performs the **mesh refinement** by:
- Adding intermediate nodes between the geometry nodes,
- Adding corresponding elements according to `elemPerBar`.

---

### 🧠 4. Element Matrices

- **`elemStiffnessMatrix.py`** — Computes **element stiffness matrices**  
- **`elemMassMatrix.py`** — Computes **element mass matrices**

➡️ Currently supports **beam elements** only.

---

### 🧭 5. Global Assembly (`globalMassStiffMatrices.py`)

Assembles the **global stiffness** and **global mass** matrices from the element contributions.

---

### 📈 6. Frequency & Mode Computation (`MECA0029_Group16_1.py`)

This file handles:
- Computation of **eigenvalues** and **eigenvectors** of the system,  
- Extraction of the **first `desiredFreqNb` modes**,  
- **Visualization** of the mode shapes.

#### ▶️ To run the modal analysis:
```bash
python ./MECA0029_Group16_1.py
```
---

### ⚖️ 7. Total Mass Verification (`MECA0029_Group16_1b.py`)

Computes the **total mass** of the structure by simulating a **rigid body translation mode** and verifying the mass distribution.

---

## 🚀 Future Improvements

- [ ] Support for **different boundary condition types** (e.g., pinned, elastic, partial constraints)
- [ ] Support for variable number of section types  
- [ ] **Refined lumped mass modeling** (non-uniform distribution, point masses)  
- [ ] Additional **element types** beyond beams  
- [ ] Local mesh refinement per beam

---



