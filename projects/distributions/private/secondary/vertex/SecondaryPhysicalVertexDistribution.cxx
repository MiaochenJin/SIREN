#include "LeptonInjector/distributions/secondary/vertex/SecondaryPhysicalVertexDistribution.h"

#include <set>                                                    // for set
#include <array>                                                  // for array
#include <cmath>                                                  // for exp
#include <tuple>                                                  // for tie
#include <string>                                                 // for bas...
#include <vector>                                                 // for vector

#include "LeptonInjector/interactions/CrossSection.h"            // for Cro...
#include "LeptonInjector/interactions/InteractionCollection.h"  // for Cro...
#include "LeptonInjector/dataclasses/InteractionRecord.h"         // for Int...
#include "LeptonInjector/dataclasses/InteractionSignature.h"      // for Int...
#include "LeptonInjector/dataclasses/Particle.h"                  // for Par...
#include "LeptonInjector/detector/DetectorModel.h"                   // for Ear...
#include "LeptonInjector/detector/Path.h"                         // for Path
#include "LeptonInjector/detector/Coordinates.h"
#include "LeptonInjector/distributions/Distributions.h"           // for Inj...
#include "LeptonInjector/geometry/Geometry.h"                     // for Geo...
#include "LeptonInjector/math/Vector3D.h"                         // for Vec...
#include "LeptonInjector/utilities/Errors.h"                      // for Sec...
#include "LeptonInjector/utilities/Random.h"                      // for LI_...

namespace LI {
namespace distributions {

using detector::DetectorPosition;
using detector::DetectorDirection;

namespace {
double log_one_minus_exp_of_negative(double x) {
    if(x < 1e-1) {
        return std::log(x) - x/2.0 + x*x/24.0 - x*x*x*x/2880.0;
    } else if(x > 3) {
        double ex = std::exp(-x);
        double ex2 = ex * ex;
        double ex3 = ex2 * ex;
        double ex4 = ex3 * ex;
        double ex5 = ex4 * ex;
        double ex6 = ex5 * ex;
        return -(ex + ex2 / 2.0 + ex3 / 3.0 + ex4 / 4.0 + ex5 / 5.0 + ex6 / 6.0);
    } else {
        return std::log(1.0 - std::exp(-x));
    }
}
}

//---------------
// class SecondaryPhysicalVertexDistribution : public VertexPositionDistribution
//---------------


void SecondaryPhysicalVertexDistribution::SampleVertex(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::SecondaryDistributionRecord & record) const {
    LI::math::Vector3D pos = record.initial_position;
    LI::math::Vector3D dir = record.direction;

    LI::math::Vector3D endcap_0 = pos;

    LI::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), std::numeric_limits<double>::infinity());
    path.ClipToOuterBounds();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(interactions->TargetTypes().begin(), interactions->TargetTypes().end());

    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = interactions->TotalDecayLength(record.record);
    LI::dataclasses::InteractionRecord fake_record = record.record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = detector_model->GetTargetMass(target);
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSectionAllFinalStates(fake_record);
        }
    }

    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);
    if(total_interaction_depth == 0) {
        throw(LI::utilities::InjectionFailure("No available interactions along path!"));
    }

    double traversed_interaction_depth;
    if(total_interaction_depth < 1e-6) {
        traversed_interaction_depth = rand->Uniform() * total_interaction_depth;
    } else {
        double exp_m_total_interaction_depth = exp(-total_interaction_depth);

        double y = rand->Uniform();
        traversed_interaction_depth = -log(y * exp_m_total_interaction_depth + (1.0 - y));
    }

    double dist = path.GetDistanceFromStartAlongPath(traversed_interaction_depth, targets, total_cross_sections, total_decay_length);
    LI::math::Vector3D vertex = path.GetFirstPoint() + dist * path.GetDirection();

    double length = (vertex - pos) * dir;
    record.SetLength(length);
}

double SecondaryPhysicalVertexDistribution::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex);

    LI::math::Vector3D endcap_0 = record.primary_initial_position;

    LI::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), std::numeric_limits<double>::infinity());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(DetectorPosition(vertex)))
        return 0.0;

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = interactions->TargetTypes();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = interactions->TotalDecayLength(record);
    LI::dataclasses::InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = detector_model->GetTargetMass(target);
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSectionAllFinalStates(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);

    path.SetPointsWithRay(path.GetFirstPoint(), path.GetDirection(), path.GetDistanceFromStartInBounds(DetectorPosition(vertex)));

    double traversed_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);

    double interaction_density = detector_model->GetInteractionDensity(path.GetIntersections(), DetectorPosition(vertex), targets, total_cross_sections, total_decay_length);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }

    return prob_density;
}

SecondaryPhysicalVertexDistribution::SecondaryPhysicalVertexDistribution() {}

std::string SecondaryPhysicalVertexDistribution::Name() const {
    return "SecondaryPhysicalVertexDistribution";
}

std::shared_ptr<SecondaryInjectionDistribution> SecondaryPhysicalVertexDistribution::clone() const {
    return std::shared_ptr<SecondaryInjectionDistribution>(new SecondaryPhysicalVertexDistribution(*this));
}

std::tuple<LI::math::Vector3D, LI::math::Vector3D> SecondaryPhysicalVertexDistribution::InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex);

    LI::math::Vector3D endcap_0 = record.primary_initial_position;

    LI::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), std::numeric_limits<double>::infinity());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(DetectorPosition(vertex)))
        return std::tuple<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));
    return std::tuple<LI::math::Vector3D, LI::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool SecondaryPhysicalVertexDistribution::equal(WeightableDistribution const & other) const {
    const SecondaryPhysicalVertexDistribution* x = dynamic_cast<const SecondaryPhysicalVertexDistribution*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool SecondaryPhysicalVertexDistribution::less(WeightableDistribution const & other) const {
    const SecondaryPhysicalVertexDistribution* x = dynamic_cast<const SecondaryPhysicalVertexDistribution*>(&other);
    return false;
}

} // namespace distributions
} // namespace LeptonInjector