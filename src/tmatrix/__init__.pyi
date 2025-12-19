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
    out_np: npt.NDArray[np.float64]
        Stores the output result.
    dim: int
        Dimension of the output array.
    mineral_property_np: npt.NDArray[np.float64]
        Contains mineral bulk modulus [Pa], shear modulus [Pa] and density [kg/m³]. Shape should be (N, 3).
    fluid_property_np: npt.NDArray[np.float64]
        Contains fluid bulk modulus [Pa] and density [kg/m³], viscosity [cP] and permeability [mD]. Shape should be (N, 4).
    phi_vector_np: npt.NDArray[np.float64]
        Porosity values array. Shape should be (N,).
    in_scenario: Literal[1, 2, 3, 4]
        1: Dual porosity, mostly rounded pores
        2: Dual porosity, little rounded pores
        3: Mixed pores
        4: Flat pores and cracks
    frequency: float
        Signal frequency [Hz]
    angle_of_sym_plane: float
        Angle of symmetry plane (0 = HTI, 90 = VTI medium) [deg]
    per_inc_con: float
        Fraction of inclusions that are connected.
    per_inc_any: float
        Fraction of inclusions that are anisotropic

    Returns
    -------
    int
        0 if success, otherwise failure.
        Result will be stored in `out_np`. Output array has shape (out_N, 4).
        Column values in order are:
            Vp: Vertical P-wave velocity [m/s]
            Vsv: Vertical polarity S-wave velocity [m/s]
            Vsh: Horizontal polarity S-wave velocity [m/s]
            Rhob [kg/m^3]
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
    """Compute TMatrix Porosity, No scenario.

    Parameters
    ----------
    out_np: npt.NDArray[np.float64]
        Stores the output result.
    out_N: int
        Dimension of the output array.
    mineral_property_np: npt.NDArray[np.float64]
        Contains mineral bulk modulus [Pa], shear modulus [Pa] and density [kg/m³]. Shape should be (N, 3).
    fluid_property_np: npt.NDArray[np.float64]
        Contains fluid bulk modulus [Pa] and density [kg/m³], viscosity [cP] and permeability [mD]. Shape should be (N, 4).
    phi_vector_np: npt.NDArray[np.float64]
        Porosity values array. Shape should be (N,) where N is the number of porosity values.
    alpha_np: npt.NDArray[np.float64]
        Aspect ratio values array. Shape should be (N,) where N is the number of aspect ratio values.
    v_np: npt.NDArray[np.float64],
        Fraction of porosity with given aspect ratio.
    alpha_size_np: npt.NDArray[np.int32],
        Number of aspect ratio values per sample,
    alpha_N: int
        Length of `alpha_np`
    frequency: float
        Signal frequency [Hz]
    angle: float
        Angle of symmetry plane (0 = HTI, 90 = VTI medium) [deg]
    inc_con_np: npt.NDArray[np.float64]
        Fraction of inclusions that are connected. Should be numpy array with a single element.
    inc_ani_np: npt.NDArray[np.float64]
        Fraction of inclusions that are anisotropic. Should be numpy array with a single element.
    inc_con_N: int
        Length of `inc_con_np` and `inc_ani_np`

    Returns
    -------
    None
        Result will be stored in `out_np`. Output array has shape (out_N, 4).
        Column values in order are:
            Vp: Vertical P-wave velocity [m/s]
            Vsv: Vertical polarity S-wave velocity [m/s]
            Vsh: Horizontal polarity S-wave velocity [m/s]
            Rhob [kg/m^3]
    """
