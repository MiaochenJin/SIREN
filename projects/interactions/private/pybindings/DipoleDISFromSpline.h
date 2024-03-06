#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/interactions/CrossSection.h"
#include "../../public/LeptonInjector/interactions/DipoleDISFromSpline.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_DipoleDISFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::interactions;

    class_<DipoleDISFromSpline, std::shared_ptr<DipoleDISFromSpline>, CrossSection> dipoledisfromspline(m, "DipoleDISFromSpline");

    dipoledisfromspline
        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::set<LI::dataclasses::ParticleType>, std::set<LI::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::vector<LI::dataclasses::ParticleType>, std::vector<LI::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, int, double, double, std::set<LI::dataclasses::ParticleType>, std::set<LI::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, std::set<LI::dataclasses::ParticleType>, std::set<LI::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, int, double, double, std::vector<LI::dataclasses::ParticleType>, std::vector<LI::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, std::vector<LI::dataclasses::ParticleType>, std::vector<LI::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(self == self)
        .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DipoleDISFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::ParticleType, double>(&DipoleDISFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DipoleDISFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::ParticleType, double, double, double, double>(&DipoleDISFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DipoleDISFromSpline::InteractionThreshold)
        .def("GetPossibleTargets",&DipoleDISFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DipoleDISFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DipoleDISFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DipoleDISFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DipoleDISFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&DipoleDISFromSpline::FinalStateProbability);
}