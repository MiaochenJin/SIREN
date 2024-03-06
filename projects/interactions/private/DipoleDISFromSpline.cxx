#include "LeptonInjector/interactions/DipoleDISFromSpline.h"

#include <map>                                             // for map, opera...
#include <set>                                             // for set, opera...
#include <array>                                           // for array
#include <cmath>                                           // for pow, log10
#include <tuple>                                           // for tie, opera...
#include <memory>                                          // for allocator
#include <string>                                          // for basic_string
#include <vector>                                          // for vector
#include <assert.h>                                        // for assert
#include <stddef.h>                                        // for size_t

#include <rk/rk.hh>                                        // for P4, Boost
#include <rk/geom3.hh>                                     // for Vector3

#include <photospline/splinetable.h>                       // for splinetable

#include "LeptonInjector/interactions/CrossSection.h"     // for CrossSection
#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/dataclasses/Particle.h"           // for Particle
#include "LeptonInjector/utilities/Random.h"               // for LI_random

namespace LI {
namespace interactions {

namespace {
///Check whether a given point in phase space is physically realizable.
///Based on equation 5 of https://arxiv.org/pdf/1707.08573.pdf
///\param x Bjorken x of the interaction
///\param y Bjorken y of the interaction
///\param E Incoming neutrino in energy in the lab frame ($E_\nu$)
///\param M Mass of the target nucleon ($M_N$)
///\param m Mass of the secondary lepton ($m_HNL$)
bool kinematicallyAllowed(double x, double y, double E, double M, double m) {
    double Q2 = 2*M*E*x*y;
    double W2 = M*M + Q2/x * (1-x);
    double Er = E*y;
    double term = M*M - W2 - 2*x*E*M - x*x*M*M + 2*Er(x*M + E,2);
    return Er*Er - W2 - term*term/(4*E*E) > 0; // equation 5
}
}

DipoleDISFromSpline::DipoleDISFromSpline() {}

DipoleDISFromSpline::DipoleDISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, double hnl_mass, std::vector<double> diople_coupling, int interaction, double target_mass, double minimum_Q2, std::set<LI::dataclasses::ParticleType> primary_types, std::set<LI::dataclasses::ParticleType> target_types, std::string units) : hnl_mass_(hnl_mass), dipole_coupling_(dipole_coupling), primary_types_(primary_types), target_types_(target_types), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    LoadFromMemory(differential_data, total_data);
    InitializeSignatures();
    SetUnits(units);
}

DipoleDISFromSpline::DipoleDISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, double hnl_mass, std::vector<double> diople_coupling, int interaction, double target_mass, double minimum_Q2, std::vector<LI::dataclasses::ParticleType> primary_types, std::vector<LI::dataclasses::ParticleType> target_types, std::string units) : hnl_mass_(hnl_mass), dipole_coupling_(dipole_coupling), primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    LoadFromMemory(differential_data, total_data);
    InitializeSignatures();
    SetUnits(units);
}

DipoleDISFromSpline::DipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> diople_coupling, int interaction, double target_mass, double minimum_Q2, std::set<LI::dataclasses::ParticleType> primary_types, std::set<LI::dataclasses::ParticleType> target_types, std::string units) : hnl_mass_(hnl_mass), dipole_coupling_(dipole_coupling), primary_types_(primary_types), target_types_(target_types), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
    SetUnits(units);
}

DipoleDISFromSpline::DipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> diople_coupling, std::set<LI::dataclasses::ParticleType> primary_types, std::set<LI::dataclasses::ParticleType> target_types, std::string units) : hnl_mass_(hnl_mass), dipole_coupling_(dipole_coupling), primary_types_(primary_types), target_types_(target_types) {
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
    SetUnits(units);
}

DipoleDISFromSpline::DipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> diople_coupling, int interaction, double target_mass, double minimum_Q2, std::vector<LI::dataclasses::ParticleType> primary_types, std::vector<LI::dataclasses::ParticleType> target_types, std::string units) : hnl_mass_(hnl_mass), dipole_coupling_(dipole_coupling), primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
    SetUnits(units);
}

DipoleDISFromSpline::DipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> diople_coupling, std::vector<LI::dataclasses::ParticleType> primary_types, std::vector<LI::dataclasses::ParticleType> target_types, std::string units) : hnl_mass_(hnl_mass), dipole_coupling_(dipole_coupling), primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()) {
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
    SetUnits(units);
}

bool DipoleDISFromSpline::equal(CrossSection const & other) const {
    const DipoleDISFromSpline* x = dynamic_cast<const DipoleDISFromSpline*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(
            target_mass_,
            hnl_mass_,
            dipole_coupling_,
            minimum_Q2_,
            signatures_,
            primary_types_,
            target_types_
            *pdf_,
            pids_)
            ==
            std::tie(
            x->target_mass_,
            x->hnl_mass_,
            x->dipole_coupling_,
            x->minimum_Q2_,
            x->signatures_,
            x->primary_types_,
            x->target_types_
            *(x->pdf_),
            x->pids_);
}

void DipoleDISFromSpline::InitializeSignatures() {
    signatures_.clear();
    for(auto primary_type : primary_types_) {
        dataclasses::InteractionSignature signature;
        signature.primary_type = primary_type;

        if(not isNeutrino(primary_type)) {
            throw std::runtime_error("This DIS implementation only supports neutrinos as primaries!");
        }

        LI::dataclasses::Particle::ParticleType neutral_lepton_product = LI::dataclasses::Particle::ParticleType::unknown;

        if(int(primary_type) > 0) {
            neutral_lepton_product = LI::dataclasses::Particle::ParticleType::NuF4;
        } 
        else {
            neutral_lepton_product = LI::dataclasses::Particle::ParticleType::NuF4Bar;
        }
        signature.secondary_types.push_back(neutral_lepton_product);

        signature.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::Hadrons);
        for(auto target_type : target_types_) {
            signature.target_type = target_type;

            signatures_.push_back(signature);

            std::pair<LI::dataclasses::Particle::ParticleType, LI::dataclasses::Particle::ParticleType> key(primary_type, target_type);
            signatures_by_parent_types_[key].push_back(signature);
        }
    }
}

double DipoleDISFromSpline::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    LI::dataclasses::ParticleType primary_type = interaction.signature.primary_type;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    double primary_energy;
    primary_energy = interaction.primary_momentum[0];
    // if we are below threshold, return 0
    if(primary_energy < InteractionThreshold(interaction))
        return 0;
    return TotalCrossSection(primary_type, primary_energy);
}

double DipoleDISFromSpline::TotalCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy) const {
    if(not primary_types_.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
    double log_energy = log10(primary_energy);

    if(log_energy < total_cross_section_.lower_extent(0)
            or log_energy > total_cross_section_.upper_extent(0)) {
        throw std::runtime_error("Interaction energy ("+ std::to_string(primary_energy) +
                ") out of cross section table range: ["
                + std::to_string(pow(10.,total_cross_section_.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10.,total_cross_section_.upper_extent(0))) + " GeV]");
    }

    int center;
    total_cross_section_.searchcenters(&log_energy, &center);
    double log_xs = total_cross_section_.ndsplineeval(&log_energy, &center, 0);
    double norm = 0;
    if (primary_type==LI::dataclasses::ParticleType::NuE || primary_type==LI::dataclasses::ParticleType::NuEBar)
        norm = std::pow(dipole_coupling_[0],2);
    else if (primary_type==LI::dataclasses::ParticleType::NuMu || primary_type==LI::dataclasses::ParticleType::NuMuBar)
        norm = std::pow(dipole_coupling_[1],2);
    else if (primary_type==LI::dataclasses::ParticleType::NuTau || primary_type==LI::dataclasses::ParticleType::NuTauBar)
        norm = std::pow(dipole_coupling_[1],2);

    return norm * unit * std::pow(10.0, log_xs);
}

// No implementation for DIS yet, just use non-target function
double DipoleDISFromSpline::TotalCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, LI::dataclasses::Particle::ParticleType target_type) const {
		return DipoleDISFromSpline::TotalCrossSection(primary_type,primary_energy);
}


double DipoleDISFromSpline::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    double primary_energy;
    primary_energy = interaction.primary_momentum[0];
    assert(interaction.signature.secondary_types.size() == 2);
    unsigned int hnl_index = (interaction.signature.secondary_types[0]==LI::dataclasses::ParticleType::NuF4 || interaction.signature.secondary_types[0]==LI::dataclasses::ParticleType::NuF4Bar) ? 0 : 1;
    unsigned int other_index = 1 - hnl_index;

    std::array<double, 4> const & mom3 = interaction.secondary_momenta[hnl_index];
    std::array<double, 4> const & mom4 = interaction.secondary_momenta[other_index];
    rk::P4 p3(geom3::Vector3(mom3[1], mom3[2], mom3[3]), interaction.secondary_masses[hnl_index]);
    rk::P4 p4(geom3::Vector3(mom4[1], mom4[2], mom4[3]), interaction.secondary_masses[other_index]);

    rk::P4 q = p1 - p3;

    double Q2 = -q.dot(q);
    double y = 1.0 - p2.dot(p3) / p2.dot(p1);
    double x = Q2 / (2.0 * p2.dot(q));

    return DifferentialCrossSection(interaction.signature.primary_type, primary_energy, x, y, Q2);
}

double DipoleDISFromSpline::DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double energy, double x, double y, double Q2) const {
    double log_energy = log10(energy);
    // check preconditions
    if(log_energy < differential_cross_section_.lower_extent(0)
            || log_energy>differential_cross_section_.upper_extent(0))
        return 0.0;
    if(x <= 0 || x >= 1)
        return 0.0;
    if(y <= 0 || y >= 1)
        return 0.0;

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    if(std::isnan(Q2)) {
        Q2 = 2.0 * energy * target_mass_ * x * y;
    }
    if(Q2 < minimum_Q2_) // cross section not calculated, assumed to be zero
        return 0;

    // cross section should be zero, but this check is missing from the original
    // CSMS calculation, so we must add it here
    if(!kinematicallyAllowed(x, y, energy, target_mass_, secondary_lepton_mass))
        return 0;

    std::array<double,3> coordinates{{log_energy, log10(x), log10(y)}};
    std::array<int,3> centers;
    if(!differential_cross_section_.searchcenters(coordinates.data(), centers.data()))
        return 0;
    double result = pow(10., differential_cross_section_.ndsplineeval(coordinates.data(), centers.data(), 0));
    assert(result >= 0);
    double norm = 0;
    if (primary_type==LI::dataclasses::ParticleType::NuE || primary_type==LI::dataclasses::ParticleType::NuEBar)
        norm = std::pow(dipole_coupling_[0],2);
    else if (primary_type==LI::dataclasses::ParticleType::NuMu || primary_type==LI::dataclasses::ParticleType::NuMuBar)
        norm = std::pow(dipole_coupling_[1],2);
    else if (primary_type==LI::dataclasses::ParticleType::NuTau || primary_type==LI::dataclasses::ParticleType::NuTauBar)
        norm = std::pow(dipole_coupling_[1],2);
    return norm * unit * result;
}

double DipoleDISFromSpline::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    // Consider implementing thershold at some point
    return 0;
}

void DipoleDISFromSpline::SampleFinalState(dataclasses::InteractionRecord& interaction, std::shared_ptr<LI::utilities::LI_random> random) const {
    // Uses Metropolis-Hastings Algorithm!
    // useful for cases where we don't know the supremum of our distribution, and the distribution is multi-dimensional
    if (differential_cross_section_.get_ndim() != 3) {
        throw std::runtime_error("I expected 3 dimensions in the cross section spline, but got " + std::to_string(differential_cross_section_.get_ndim()) +". Maybe your fits file doesn't have the right 'INTERACTION' key?");
    }

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * tinteraction.secondary_momentarget_mass_ + 2 * target_mass_ * primary_energy;
    // double s = std::pow(rk::invMass(p1, p2), 2);

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    p1_lab = p1;
    p2_lab = p2;
    primary_energy = p1_lab.e();

    unsigned int hnl_index = (interaction.signature.secondary_types[0]==LI::dataclasses::ParticleType::NuF4 || interaction.signature.secondary_types[0]==LI::dataclasses::ParticleType::NuF4Bar) ? 0 : 1;
    unsigned int other_index = 1 - hnl_index;
    double m = hnl_mass_;

    double m1 = interaction.primary_mass;
    double m3 = m;
    double E1_lab = p1_lab.e();
    double E2_lab = p2_lab.e();


    // The out-going particle always gets at least enough energy for its rest mass
    double yMax = 1 - m / primary_energy;
    double logYMax = log10(yMax);

    // The minimum allowed value of y occurs when x = 1 and Q is minimized
    double yMin = minimum_Q2_ / (2 * E1_lab * E2_lab);
    double logYMin = log10(yMin);
    // The minimum allowed value of x occurs when y = yMax and Q is minimized
    // double xMin = minimum_Q2_ / ((s - target_mass_ * target_mass_) * yMax);
    double xMin = minimum_Q2_ / (2 * E1_lab * E2_lab * yMax);
    double logXMin = log10(xMin);

    bool accept;

    // kin_vars and its twin are 3-vectors containing [nu-energy, Bjorken X, Bjorken Y]
    std::array<double,3> kin_vars, test_kin_vars;

    // centers of the cross section spline tales.
    std::array<int,3> spline_table_center, test_spline_table_center;

    // values of cross_section from the splines.  By * Bx * Spline(E,x,y)
    double cross_section, test_cross_section;

    // No matter what, we're evaluating at this specific energy.
    kin_vars[0] = test_kin_vars[0] = log10(primary_energy);

    // check preconditions
    if(kin_vars[0] < differential_cross_section_.lower_extent(0)
            || kin_vars[0] > differential_cross_section_.upper_extent(0))
        throw std::runtime_error("Interaction energy out of cross section table range: ["
                + std::to_string(pow(10.,differential_cross_section_.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10.,differential_cross_section_.upper_extent(0))) + " GeV]");

    // sample an intial point
    do {
        // rejection sample a point which is kinematically allowed by calculation limits
        double trialQ;
        do {
            kin_vars[1] = random->Uniform(logXMin,0);
            kin_vars[2] = random->Uniform(logYMin,logYMax);
            trialQ = (2 * E1_lab * E2_lab) * pow(10., kin_vars[1] + kin_vars[2]);
        } while(trialQ<minimum_Q2_ || !kinematicallyAllowed(pow(10., kin_vars[1]), pow(10., kin_vars[2]), primary_energy, target_mass_, m));

        accept = true;
        //sanity check: demand that the sampled point be within the table extents
        if(kin_vars[1] < differential_cross_section_.lower_extent(1)
                || kin_vars[1] > differential_cross_section_.upper_extent(1)) {
            accept = false;
        }
        if(kin_vars[2] < differential_cross_section_.lower_extent(2)
                || kin_vars[2] > differential_cross_section_.upper_extent(2)) {
            accept = false;
        }

        if(accept) {
            // finds the centers in the cross section spline table, returns true if it's successful
            // also sets the centers
            accept = differential_cross_section_.searchcenters(kin_vars.data(),spline_table_center.data());
        }
    } while(!accept);

    //TODO: better proposal distribution?
    double measure = pow(10., kin_vars[1] + kin_vars[2]); // Bx * By

    // Bx * By * xs(E, x, y)
    // evalutates the differential spline at that point
    cross_section = measure*pow(10., differential_cross_section_.ndsplineeval(kin_vars.data(), spline_table_center.data(), 0));

    // this is the magic part. Metropolis Hastings Algorithm.
    // MCMC method!
    const size_t burnin = 40; // converges to the correct distribution over multiple samplings.
    // big number means more accurate, but slower
    for(size_t j = 0; j <= burnin; j++) {
        // repeat the sampling from above to get a new valid point
        double trialQ;
        do {
            test_kin_vars[1] = random->Uniform(logXMin, 0);
            test_kin_vars[2] = random->Uniform(logYMin, logYMax);
            trialQ = (2 * E1_lab * E2_lab) * pow(10., test_kin_vars[1] + test_kin_vars[2]);
        } while(trialQ < minimum_Q2_ || !kinematicallyAllowed(pow(10., test_kin_vars[1]), pow(10., test_kin_vars[2]), primary_energy, target_mass_, m));

        accept = true;
        if(test_kin_vars[1] < differential_cross_section_.lower_extent(1)
                || test_kin_vars[1] > differential_cross_section_.upper_extent(1))
            accept = false;
        if(test_kin_vars[2] < differential_cross_section_.lower_extent(2)
                || test_kin_vars[2] > differential_cross_section_.upper_extent(2))
            accept = false;
        if(!accept)
            continue;

        accept = differential_cross_section_.searchcenters(test_kin_vars.data(), test_spline_table_center.data());
        if(!accept)
            continue;

        double measure = pow(10., test_kin_vars[1] + test_kin_vars[2]);
        double eval = differential_cross_section_.ndsplineeval(test_kin_vars.data(), test_spline_table_center.data(), 0);
        if(std::isnan(eval))
            continue;
        test_cross_section = measure * pow(10., eval);

        double odds = (test_cross_section / cross_section);
        accept = (cross_section == 0 || (odds > 1.) || random->Uniform(0, 1) < odds);

        if(accept) {
            kin_vars = test_kin_vars;
            cross_section = test_cross_section;
        }
    }
    double final_x = pow(10., kin_vars[1]);
    double final_y = pow(10., kin_vars[2]);

    interaction.interaction_parameters.clear();
    record.interaction_parameters["energy"] = E1_lab;
    record.interaction_parameters["bjorken_x"] = final_x;
    record.interaction_parameters["bjorken_y"] = final_y;

    double Q2 = 2 * E1_lab * E2_lab * pow(10.0, kin_vars[1] + kin_vars[2]);
    double p1x_lab = std::sqrt(p1_lab.px() * p1_lab.px() + p1_lab.py() * p1_lab.py() + p1_lab.pz() * p1_lab.pz());
    double pqx_lab = (m1*m1 + m3*m3 + 2 * p1x_lab * p1x_lab + Q2 + 2 * E1_lab * E1_lab * (final_y - 1)) / (2.0 * p1x_lab);
    double momq_lab = std::sqrt(m1*m1 + p1x_lab*p1x_lab + Q2 + E1_lab * E1_lab * (final_y * final_y - 1));
    double pqy_lab = std::sqrt(momq_lab*momq_lab - pqx_lab *pqx_lab);
    double Eq_lab = E1_lab * final_y;

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 p1_mom = p1_lab.momentum();
    geom3::UnitVector3 p1_lab_dir = p1_mom.direction();
    geom3::Rotation3 x_to_p1_lab_rot = geom3::rotationBetween(x_dir, p1_lab_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(p1_lab_dir, phi);

    rk::P4 pq_lab(Eq_lab, geom3::Vector3(pqx_lab, pqy_lab, 0));
    pq_lab.rotate(x_to_p1_lab_rot);
    pq_lab.rotate(rand_rot);

    rk::P4 p3_lab((p1_lab - pq_lab).momentum(), m3);
    rk::P4 p4_lab = p2_lab + pq_lab;

    rk::P4 p3;
    rk::P4 p4;
    p3 = p3_lab;
    p4 = p4_lab;

    std::vector<LI::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    LI::dataclasses::SecondaryParticleRecord & hnl = secondaries[hnl_index];
    LI::dataclasses::SecondaryParticleRecord & other = secondaries[other_index];


    hnl.SetFourMomentum({p3.e(), p3.px(), p3.py(), p3.pz()});
    hnl.SetMass(p3.m());
    hnl.SetHelicity(-1.0 * record.primary_helicity); // assume helicity flipping

    other.SetFourMomentum({p4.e(), p4.px(), p4.py(), p4.pz()});
    other.SetMass(p4.m());
    other.SetHelicity(record.target_helicity);
}

double DipoleDISFromSpline::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<LI::dataclasses::Particle::ParticleType> DipoleDISFromSpline::GetPossiblePrimaries() const {
    return std::vector<LI::dataclasses::Particle::ParticleType>(primary_types_.begin(), primary_types_.end());
}

std::vector<LI::dataclasses::Particle::ParticleType> DipoleDISFromSpline::GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const {
    return std::vector<LI::dataclasses::Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> DipoleDISFromSpline::GetPossibleSignatures() const {
    return std::vector<dataclasses::InteractionSignature>(signatures_.begin(), signatures_.end());
}

std::vector<LI::dataclasses::Particle::ParticleType> DipoleDISFromSpline::GetPossibleTargets() const {
    return std::vector<LI::dataclasses::Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> DipoleDISFromSpline::GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const {
    std::pair<LI::dataclasses::Particle::ParticleType, LI::dataclasses::Particle::ParticleType> key(primary_type, target_type);
    if(signatures_by_parent_types_.find(key) != signatures_by_parent_types_.end()) {
        return signatures_by_parent_types_.at(key);
    } else {
        return std::vector<dataclasses::InteractionSignature>();
    }
}

std::vector<std::string> DipoleDISFromSpline::DensityVariables() const {
    return std::vector<std::string>{"Bjorken x", "Bjorken y"};
}

} // namespace interactions
} // namespace LI
