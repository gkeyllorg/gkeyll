# Surrogate Models Overview

### `nn_model_sheath_bc_conv_MPE.kann`

A neural network model trained on simulation data from **GYRAZE**. It predicts the electron cutting velocity profile based on local plasma and magnetic parameters at the magnetic presheath entrance (MPE).

#### Inputs

* **`alpha`**: Angle between the magnetic field and the wall tangent.
* **`gamma`**: Ratio between the electron Larmor radius and the Debye length:

$$\gamma = \frac{1}{B} \sqrt{\frac{n_e}{\varepsilon_0}}$$


* **`phi`**: Normalized potential at the wall:

$$\phi_{\text{norm}} = \frac{e \phi}{T_e}$$



#### Outputs

The model outputs a vector of **40 values**:

* **Indices 0–19**: Normalized electron cutting velocities $\hat v_i$ [$v_\text{th}$]
* **Indices 20–39**: Corresponding normalized magnetic moments $\hat \mu_i$ [$T_e/ B$]

#### References & Documentation

For full details on the architecture, training process, and C interface, visit the repository:
[Antoinehoff/gyraze_surrogate_c_interface](https://github.com/Antoinehoff/gyraze_surrogate_c_interface)