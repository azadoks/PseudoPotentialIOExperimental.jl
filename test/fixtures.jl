import PseudoPotentialIOExperimental: resolve_family

TEST_FAMILIES = (
    "hgh_lda_hgh",
    "hgh_pbe_upf",
    "gbrv_lda_1.5_upf",
    "sg15_2022.02.06_upf",
    "sssp_pbe_precision_1.1.2_upf",
    "pd_nc_sr_pbe_standard_0.4.1_upf",
    "pd_nc_sr_pbe_standard_0.4.1_psp8",
    "pd_nc_fr_pbesol_standard_0.4_psp8",
)
TEST_FILEPATHS = []
for family in TEST_FAMILIES
    family_dir = resolve_family(family)
    psp_filenames = list_family_psps(family)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(TEST_FILEPATHS, psp_filepaths)
end

UPF_FAMILIES = (
    "gbrv_pbe_1.5_upf",
    "hgh_lda_upf",
    "sg15_2022.02.06_upf",
    "sssp_pbe_efficiency_1.1.2_upf",
)
UPF_FILEPATHS = []
for family in UPF_FAMILIES
    family_dir = resolve_family(family)
    psp_filenames = list_family_psps(family)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(UPF_FILEPATHS, psp_filepaths)
end

HGH_FAMILIES = (
    "hgh_lda_hgh",
    "hgh_pbe_hgh",
)
HGH_FILEPATHS = []
for family in HGH_FAMILIES
    family_dir = resolve_family(family)
    psp_filenames = list_family_psps(family)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(HGH_FILEPATHS, psp_filepaths)
end

PSP8_FAMILIES = (
    "pd_nc_sr_pbe_standard_0.4.1_psp8",
    "pd_nc_fr_pbesol_standard_0.4_psp8",
)
PSP8_FILEPATHS = []
for family in PSP8_FAMILIES
    family_dir = resolve_family(family)
    psp_filenames = list_family_psps(family)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(PSP8_FILEPATHS, psp_filepaths)
end

UPF1_CASES = [
    ("gbrv_lda_1.5_upf", "ag_lda_v1.4.uspp.F.UPF"),
    ("sssp_pbe_precision_1.1.2_upf", "B_pbe_v1.01.uspp.F.UPF"),
    ("gbrv_pbesol_1.5_upf", "si_pbesol_v1.uspp.F.UPF"),
    ("gbrv_pbe_1.5_upf", "mg_pbe_v1.4.uspp.F.UPF"),
]
UPF1_CASE_FILEPATHS = Dict()
for (family, filename) in UPF1_CASES
    if !haskey(UPF1_CASE_FILEPATHS, filename)
        UPF1_CASE_FILEPATHS[filename] = joinpath(resolve_family(family), filename)
    else
        error("UPF1_CASE_FILEPATHS already contains a case for $filename")
    end
end

UPF2_CASES = [
    ("pd_nc_fr_pbesol_stringent_0.4_upf", "Mg.upf"),
    ("sssp_pbe_efficiency_1.1.2_upf", "Si.pbe-n-rrkjus_psl.1.0.0.UPF"),
    ("sssp_pbe_efficiency_1.1.2_upf", "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ("sssp_pbe_efficiency_1.1.2_upf", "Dy.GGA-PBE-paw-v1.0.UPF"),
    ("hgh_pbe_upf", "He.pbe-hgh.UPF"),
    ("hgh_pbe_upf", "Al.pbe-hgh.UPF"),
    ("pd_nc_sr_pbe_standard_0.4.1_upf", "H.upf"),
    ("pd_nc_sr_pbe_standard_0.4.1_upf", "Fe.upf"),
]
UPF2_CASE_FILEPATHS = Dict()
for (family, filename) in UPF2_CASES
    if !haskey(UPF2_CASE_FILEPATHS, filename)
        UPF2_CASE_FILEPATHS[filename] = joinpath(resolve_family(family), filename)
    else
        error("UPF2_CASE_FILEPATHS already contains a case for $filename")
    end
end

PSP8_CASES = [
     ("pd_nc_sr_pbe_standard_0.4.1_psp8", "H.psp8"),
     ("pd_nc_sr_pbe_standard_0.4.1_psp8", "Fe.psp8"),
     ("pd_nc_fr_pbesol_standard_0.4_psp8", "Ti.psp8"),
     ("pd_nc_sr_pbesol_stringent_0.4.1_psp8", "Zn.psp8"),
]
PSP8_CASE_FILEPATHS = Dict()
for (family, filename) in PSP8_CASES
    if !haskey(PSP8_CASE_FILEPATHS, filename)
        PSP8_CASE_FILEPATHS[filename] = joinpath(resolve_family(family), filename)
    else
        error("PSP8_CASE_FILEPATHS already contains a case for $filename")
    end
end

HGH_CASES = [
    ("hgh_lda_hgh", "c-q4.hgh"),
    ("hgh_lda_hgh", "ni-q18.hgh"),
]
HGH_CASE_FILEPATHS = Dict()
for (family, filename) in HGH_CASES
    if !haskey(HGH_CASE_FILEPATHS, filename)
        HGH_CASE_FILEPATHS[filename] = joinpath(resolve_family(family), filename)
    else
        error("HGH_CASE_FILEPATHS already contains a case for $filename")
    end
end

UPF2_PSP8_ELEMENTS = ("Ag", "Al", "Ar", "As", "Au", "B", "Ba", "Be", "Bi", "Br", "C", "Ca",
                      # "Cr", "Cu", # TODO Chromium and Copper agreement is terrible for some reason
                      "Cd", "Cl", "Co", "Cs", "F", "Fe", "Ga", "Ge", "H", "He",
                      "Hf", "Hg", "I", "In", "Ir", "K", "Kr", "La", "Li", "Lu", "Mg", "Mn",
                      "Mo", "N", "Na", "Nb", "Ne", "Ni", "O", "Os", "P", "Pb", "Pd", "Po",
                      "Pt", "Rb", "Re", "Rh", "Rn", "Ru", "S", "Sb", "Sc", "Se", "Si", "Sn",
                      "Sr", "Ta", "Tc", "Te", "Ti", "Tl", "V", "W", "Xe", "Y", "Zn", "Zr")
UPF2_PSP8_FILEPATHS = map(UPF2_PSP8_ELEMENTS) do element
    (
        joinpath(resolve_family("pd_nc_sr_pbe_standard_0.4.1_upf"), "$(element).upf"),
        joinpath(resolve_family("pd_nc_sr_pbe_standard_0.4.1_psp8"), "$(element).psp8")
    )
end
