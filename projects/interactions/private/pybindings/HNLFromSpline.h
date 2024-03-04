#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/interactions/CrossSection.h"
#include "../../public/LeptonInjector/interactions/HNLFromSpline.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_HNLFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::interactions;

    class_<HNLFromSpline, std::shared_ptr<HNLFromSpline>, CrossSection> disfromspline(m, "HNLFromSpline");

    disfromspline
        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::set<LI::dataclasses::ParticleType>, std::set<LI::dataclasses::ParticleType>>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::vector<LI::dataclasses::ParticleType>, std::vector<LI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::set<LI::dataclasses::ParticleType>, std::set<LI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, std::set<LI::dataclasses::ParticleType>, std::set<LI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::vector<LI::dataclasses::ParticleType>, std::vector<LI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, std::vector<LI::dataclasses::ParticleType>, std::vector<LI::dataclasses::ParticleType>>())
        .def(self == self)
        .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&HNLFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::ParticleType, double>(&HNLFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::ParticleType, double, LI::dataclasses::ParticleType>(&HNLFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&HNLFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<double, double, double, double>(&HNLFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&HNLFromSpline::InteractionThreshold)
        .def("GetPossibleTargets",&HNLFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&HNLFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&HNLFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&HNLFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&HNLFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&HNLFromSpline::FinalStateProbability);
}

