
#include <vector>
#include <set>

#include "../../public/LeptonInjector/dataclasses/Particle.h"
#include "../../public/LeptonInjector/dataclasses/ParticleID.h"
#include "../../public/LeptonInjector/dataclasses/ParticleType.h"
#include "../../public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../public/LeptonInjector/dataclasses/InteractionTree.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace pybind11;

PYBIND11_MODULE(dataclasses,m) {
  using namespace LI::dataclasses;

  class_<Particle, std::shared_ptr<Particle>> particle(m, "Particle");

  particle.def(init<>())
          .def(init<Particle const &>())
          .def(init<ParticleID, ParticleType, double, std::array<double, 4>, std::array<double, 3>, double, double>())
          .def(init<ParticleType, double, std::array<double, 4>, std::array<double, 3>, double, double>())
          .def_readwrite("id",&Particle::id)
          .def_readwrite("type",&Particle::type)
          .def_readwrite("mass",&Particle::mass)
          .def_readwrite("momentum",&Particle::momentum)
          .def_readwrite("position",&Particle::position)
          .def_readwrite("length",&Particle::length)
          .def_readwrite("helicity",&Particle::helicity)
          .def("GenerateID",&Particle::GenerateID);

  enum_<ParticleType>(particle, "ParticleType")
#define X(a, b) .value( #a , ParticleType:: a )
#include "../../public/LeptonInjector/dataclasses/ParticleTypes.def"
#undef X
          .export_values();

  class_<InteractionSignature, std::shared_ptr<InteractionSignature>>(m, "InteractionSignature")
          .def(init<>())
          .def_readwrite("primary_type",&InteractionSignature::primary_type)
          .def_readwrite("target_type",&InteractionSignature::target_type)
          .def_readwrite("secondary_types",&InteractionSignature::secondary_types);

    class_<PrimaryDistributionRecord, std::shared_ptr<PrimaryDistributionRecord>>(m, "PrimaryDistributionRecord")
        .def(init<InteractionRecord const &>())
        .def_property_readonly("record",
            [](LI::dataclasses::PrimaryDistributionRecord const & pdr) {LI::dataclasses::InteractionRecord ir = pdr.record; return ir;})
        .def_property_readonly("signature",
            [](LI::dataclasses::PrimaryDistributionRecord const & pdr) {LI::dataclasses::InteractionSignature is = pdr.signature; return is;})
        .def_property_readonly("id",
            [](LI::dataclasses::PrimaryDistributionRecord const & pdr) {LI::dataclasses::ParticleID id = pdr.id; return id;})
        .def_property_readonly("type",
            [](LI::dataclasses::PrimaryDistributionRecord const & pdr) {LI::dataclasses::ParticleType pt = pdr.type; return pt;})
        .def("GetParticle", &PrimaryDistributionRecord::GetParticle)
        .def("SetParticle", &PrimaryDistributionRecord::SetParticle)
        .def_property("mass", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetMass)), &PrimaryDistributionRecord::SetMass)
        .def_property("energy", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetEnergy)), &PrimaryDistributionRecord::SetEnergy)
        .def_property("kinetic_energy", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetKineticEnergy)), &PrimaryDistributionRecord::SetKineticEnergy)
        .def_property("direction", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetDirection)), &PrimaryDistributionRecord::SetDirection)
        .def_property("three_momentum", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetThreeMomentum)), &PrimaryDistributionRecord::SetThreeMomentum)
        .def_property("four_momentum", ((std::array<double, 4> (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetFourMomentum)), &PrimaryDistributionRecord::SetFourMomentum)
        .def_property("length", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetLength)), &PrimaryDistributionRecord::SetLength)
        .def_property("initial_position", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetInitialPosition)), &PrimaryDistributionRecord::SetInitialPosition)
        .def_property("interaction_vertex", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetInteractionVertex)), &PrimaryDistributionRecord::SetInteractionVertex)
        .def_property("helicity", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetHelicity)), &PrimaryDistributionRecord::SetHelicity)
        .def("Finalize", &PrimaryDistributionRecord::Finalize);

    class_<SecondaryParticleRecord, std::shared_ptr<SecondaryParticleRecord>>(m, "SecondaryParticleRecord")
        .def(init<InteractionRecord const &, size_t>())
        .def_property_readonly("id",
            [](LI::dataclasses::SecondaryParticleRecord const & spr) {LI::dataclasses::ParticleID id = spr.id; return id;})
        .def_property_readonly("type",
            [](LI::dataclasses::SecondaryParticleRecord const & spr) {LI::dataclasses::ParticleType pt = spr.type; return pt;})
        .def_property_readonly("initial_position",
            [](LI::dataclasses::SecondaryParticleRecord const & spr) {std::array<double, 3> ip = spr.initial_position; return ip;})
        .def("GetParticle", &SecondaryParticleRecord::GetParticle)
        .def("SetParticle", &SecondaryParticleRecord::SetParticle)
        .def_property("mass", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetMass)), &SecondaryParticleRecord::SetMass)
        .def_property("energy", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetEnergy)), &SecondaryParticleRecord::SetEnergy)
        .def_property("kinetic_energy", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetKineticEnergy)), &SecondaryParticleRecord::SetKineticEnergy)
        .def_property("direction", ((std::array<double, 3> const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetDirection)), &SecondaryParticleRecord::SetDirection)
        .def_property("three_momentum", ((std::array<double, 3> const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetThreeMomentum)), &SecondaryParticleRecord::SetThreeMomentum)
        .def_property("four_momentum", ((std::array<double, 4> (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetFourMomentum)), &SecondaryParticleRecord::SetFourMomentum)
        .def_property("helicity", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetHelicity)), &SecondaryParticleRecord::SetHelicity)
        .def("Finalize", &SecondaryParticleRecord::Finalize);

  class_<InteractionRecord, std::shared_ptr<InteractionRecord>>(m, "InteractionRecord")
          .def(init<>())
          .def_readwrite("signature",&InteractionRecord::signature)
          .def_readwrite("primary_mass",&InteractionRecord::primary_mass)
          .def_readwrite("primary_momentum",&InteractionRecord::primary_momentum)
          .def_readwrite("primary_helicity",&InteractionRecord::primary_helicity)
          .def_readwrite("target_mass",&InteractionRecord::target_mass)
          .def_readwrite("target_helicity",&InteractionRecord::target_helicity)
          .def_readwrite("interaction_vertex",&InteractionRecord::interaction_vertex)
          .def_readwrite("secondary_masses",&InteractionRecord::secondary_masses)
          .def_readwrite("secondary_momenta",&InteractionRecord::secondary_momenta)
          .def_readwrite("secondary_helicities",&InteractionRecord::secondary_helicities)
          .def_readwrite("interaction_parameters",&InteractionRecord::interaction_parameters);

  class_<InteractionTreeDatum, std::shared_ptr<InteractionTreeDatum>>(m, "InteractionTreeDatum")
          .def(init<InteractionRecord&>())
          .def_readwrite("record",&InteractionTreeDatum::record)
          .def_readwrite("parent",&InteractionTreeDatum::parent)
          .def_readwrite("daughters",&InteractionTreeDatum::daughters)
          .def("depth",&InteractionTreeDatum::depth);

  class_<InteractionTree, std::shared_ptr<InteractionTree>>(m, "InteractionTree")
          .def(init<>())
          .def_readwrite("tree",&InteractionTree::tree)
          .def("add_entry",static_cast<std::shared_ptr<InteractionTreeDatum> (InteractionTree::*)(InteractionTreeDatum&,std::shared_ptr<InteractionTreeDatum>)>(&InteractionTree::add_entry))
          .def("add_entry",static_cast<std::shared_ptr<InteractionTreeDatum> (InteractionTree::*)(InteractionRecord&,std::shared_ptr<InteractionTreeDatum>)>(&InteractionTree::add_entry));

}
