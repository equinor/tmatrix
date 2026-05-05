from typing import Literal

import numpy as np
import numpy.typing as npt

def tmatrix_porosity(
    out_np: npt.NDArray[np.float64],
    dim: int,
    mineral_property_np: npt.NDArray[np.float64],
    fluid_property_np: npt.NDArray[np.float64],
    phi_vector_np: npt.NDArray[np.float64],
    in_scenario: Literal[1, 2, 3, 4],
    frequency: float,
    angle_of_sym_plane: float,
    per_inc_con: float,
    per_inc_any: float,
) -> int:
    """Compute TMatrix Porosity.

    Parameters
    ----------
    out_np
        Stores the output result. Shape should be (dim, 4).
    dim
        Dimension of the output array.
    mineral_property_np
        Contains mineral bulk modulus [Pa], shear modulus [Pa] and density
        [kg/m³]. Shape should be (dim, 3).
    fluid_property_np
        Contains fluid bulk modulus [Pa] and density [kg/m³], viscosity [cP]
        and permeability [mD]. Shape should be (dim, 4).
    phi_vector_np
        Porosity values array. Shape should be (dim,).
    in_scenario
        Pore-shape scenario:

        - 1: Dual porosity, mostly rounded pores
        - 2: Dual porosity, little rounded pores
        - 3: Mixed pores
        - 4: Flat pores and cracks
    frequency
        Signal frequency [Hz].
    angle_of_sym_plane
        Angle of symmetry plane (0 = HTI, 90 = VTI medium) [deg].
    per_inc_con
        Fraction of inclusions that are connected.
    per_inc_any
        Fraction of inclusions that are anisotropic.

    Returns
    -------
    Always ``0``. The result is stored in ``out_np`` with shape
    (dim, 4). Columns in order are:

    - Vp: Vertical P-wave velocity [m/s]
    - Vsh: Horizontal polarity S-wave velocity [m/s]
    - Vsv: Vertical polarity S-wave velocity [m/s]
    - Rhob: Effective density [kg/m³]
    """

def tmatrix_porosity_noscenario(
    out_np: npt.NDArray[np.float64],
    out_N: int,
    mineral_property_np: npt.NDArray[np.float64],
    fluid_property_np: npt.NDArray[np.float64],
    phi_vector_np: npt.NDArray[np.float64],
    alpha_np: npt.NDArray[np.float64],
    v_np: npt.NDArray[np.float64],
    alpha_size_np: npt.NDArray[np.int32],
    alpha_N: int,
    frequency: float,
    angle: float,
    inc_con_np: npt.NDArray[np.float64],
    inc_ani_np: npt.NDArray[np.float64],
    inc_con_N: int,
) -> None:
    """Compute TMatrix Porosity, no scenario.

    Parameters
    ----------
    out_np
        Stores the output result. Shape should be (out_N, 4).
    out_N
        Dimension of the output array.
    mineral_property_np
        Contains mineral bulk modulus [Pa], shear modulus [Pa] and density
        [kg/m³]. Shape should be (out_N, 3).
    fluid_property_np
        Contains fluid bulk modulus [Pa] and density [kg/m³], viscosity [cP]
        and permeability [mD]. Shape should be (out_N, 4).
    phi_vector_np
        Porosity values array. Shape should be (out_N,).
    alpha_np
        Aspect ratio values array. Shape should be (N,) where N is the number
        of aspect ratio values.
    v_np
        Fraction of porosity with given aspect ratio.
    alpha_size_np
        Number of aspect ratio values per sample.
    alpha_N
        Length of ``alpha_np``.
    frequency
        Signal frequency [Hz].
    angle
        Angle of symmetry plane (0 = HTI, 90 = VTI medium) [deg].
    inc_con_np
        Fraction of inclusions that are connected. Should be a numpy array
        with a single element.
    inc_ani_np
        Fraction of inclusions that are anisotropic. Should be a numpy array
        with a single element.
    inc_con_N
        Length of ``inc_con_np`` and ``inc_ani_np``.

    Returns
    -------
    The result is stored in ``out_np`` with shape (out_N, 4). Columns in
    order are:

    - Vp: Vertical P-wave velocity [m/s]
    - Vsh: Horizontal polarity S-wave velocity [m/s]
    - Vsv: Vertical polarity S-wave velocity [m/s]
    - Rhob: Effective density [kg/m³]
    """
