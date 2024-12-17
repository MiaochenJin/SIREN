import os
import siren
from siren._util import GenerateEvents,SaveEvents
import nuflux
import numpy as np

# Number of events to inject
events_to_inject = int(200)

# Experiment to run
experiment = "IceCube"
detector_model = siren.utilities.load_detector(experiment)

# some environmental variables
date = '1216'
gen_emin = '1e5'
gen_emax = '1e7'
nu_type = "NuE"
current_type = 'cc'
physical_flux = "astro" # this option needs tweaking
injection = 'ranged'

# current type string
if current_type == "cc":
    current_str = "_CC"
    int_type = 1
else:
    current_str = "_NC"
    int_type = 2

# neutrino type string
if nu_type == "NuE":
    primary_type = siren.dataclasses.Particle.ParticleType.NuE
elif nu_type == "NuMu":
    primary_type = siren.dataclasses.Particle.ParticleType.NuMu

expname = "{}_{}{}_{}_nuclear_{}_{}-{}".format(date, nu_type, current_str, physical_flux, injection, gen_emin, gen_emax)

# Cross-section model to use
hydrogen_model = "PDF4LHC21_charm"
oxygen_model = "EPPS21_charm"

# Load the cross-section model
hydrogen_process, _ = siren.utilities.load_processes(
    hydrogen_model,
    primary_types=[primary_type],
    target_types=[siren.dataclasses.Particle.ParticleType.HNucleus],
    isoscalar=False,
    process_types=[current_type]
)
oxygen_process, _ = siren.utilities.load_processes(
    oxygen_model,
    primary_types=[primary_type],
    target_types=[siren.dataclasses.Particle.ParticleType.O16Nucleus],
    isoscalar=False,
    process_types=[current_type]
)
# combine the two dictionaries
primary_processes = {key: [hydrogen_process[key][0], oxygen_process[key][0]] for key in hydrogen_process}
# Extract the primary cross-sections for the primary type
primary_cross_sections = primary_processes[primary_type]

# Set up the Injector
injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_cross_sections

# position distribution for volume injection
if injection == "volume":
    # position_distribution = controller.GetCylinderVolumePositionDistributionFromSector("icecube")
    # primary_injection_distributions["position"] = position_distribution
    muon_range_func = siren.distributions.LeptonDepthFunction()
    position_distribution = siren.distributions.ColumnDepthPositionDistribution(
        600, 600.0, muon_range_func, set(controller.GetDetectorModelTargets()[0])
    )
elif injection == "ranged":
    # position distribution for ranged injection
    muon_range_func = siren.distributions.LeptonDepthFunction()
    position_distribution = siren.distributions.ColumnDepthPositionDistribution(
        600, 600.0, muon_range_func, set(controller.GetDetectorModelTargets()[0])
    )

# Directly set the distributions
injector.primary_injection_distributions = [
    siren.distributions.PrimaryMass(0),  # Mass distribution
    siren.distributions.PowerLaw(2, float(gen_emin), float(gen_emax)),  # Energy distribution
    siren.distributions.IsotropicDirection(),  # Direction distribution
    position_distribution
]

# Generate events
events,gen_times = GenerateEvents(injector)

# Set up the Weighter for event weighting (without position distribution)
weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_cross_sections
if physical_flux == 'astro':
    physical_dist = siren.distributions.PowerLaw(2.58, 1e3, 1e7)
elif physical_flux == 'atmos':
    # make an atmospheric flux
    flux = nuflux.makeFlux('honda2006')
    erange = np.logspace(2,6,100)
    erange_atmo = np.logspace(2,6,100)
    cosrange = np.linspace(0,1,100)
    atmo_flux_tables = {}
    if nu_type == "NuMu":
        particle = nuflux.NuMu
    elif nu_type == 'NuE':
        particle = nuflux.NuE
    siren_key = siren.dataclasses.Particle.ParticleType(int(particle))
    atmo_flux_tables[siren_key] = np.zeros(len(erange))
    for i,e in enumerate(erange):
        f = flux.getFlux(particle,e,cosrange)
        atmo_flux_tables[siren_key][i] += 0.01*np.sum(f) * 1e4 * 2 * np.pi
    physical_dist = siren.distributions.TabulatedFluxDistribution(erange_atmo,atmo_flux_tables[primary_type],True)
    
weighter.primary_physical_distributions = [
    physical_dist,  # Energy distribution
    siren.distributions.IsotropicDirection()  # Direction distribution
]

outdir = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/miaochenjin/DBSearch/SIREN_outputs/"
savedir = os.path.join(outdir, expname)
SaveEvents(events,weighter,gen_times,output_filename="{}/{}_".format(savedir, expname))
