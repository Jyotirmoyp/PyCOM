# PyCOM

PyCOM is a Python script designed for **post-processing spherical CitcomS simulation results**.
It enables visualization and analysis of model outputs such as scalar fields and velocity vectors.

---

### Requirements

* Python 3.x
* Required Python libraries (depending on your implementation), for example:

  * `numpy`
  * `PyGMT`


Make sure all required dependencies are installed before running the script.

---

### Configuration

Before running the script, you must configure the following parameters **at the beginning of the Python file**:

* **Model name** – the name of the CitcomS model directory.
* **Number of processors** – the number of processors used during the CitcomS simulation.

Example:

```python
model_name = "example_model"
nproc = 64
```

These parameters allow the script to locate and correctly read the distributed output files.

---

### Plotting

To generate plots at a specific depth, you must specify the **layer index** corresponding to the desired depth level in the model.

Example:

```python
layer_index = 10
```

The layer index determines which radial layer of the spherical model will be visualized.

---

### Visualization Settings

Visualization parameters can be adjusted depending on your analysis needs:

* **Color maps** for scalar fields (e.g., temperature, viscosity)
* **Velocity vectors** (density, scaling, arrow style)
* **Plot limits and scaling**
* **Figure size and resolution**
* **Text modifications**

These settings can be modified directly in the plotting section of the script.


---

### Notes

* Ensure the CitcomS output files are located in the expected directory structure.


