
#include <vector>

#include "../../public/LeptonInjector/interactions/CrossSection.h"
#include "../../public/LeptonInjector/interactions/NeutrissimoDecay.h"
#include "../../public/LeptonInjector/interactions/InteractionCollection.h"
#include "../../public/LeptonInjector/interactions/DISFromSpline.h"
#include "../../public/LeptonInjector/interactions/DipoleDISFromSpline.h"
#include "../../public/LeptonInjector/interactions/HNLFromSpline.h"
#include "../../public/LeptonInjector/interactions/DipoleFromTable.h"
#include "../../public/LeptonInjector/interactions/DarkNewsCrossSection.h"
#include "../../public/LeptonInjector/interactions/DarkNewsDecay.h"

#include "./CrossSection.h"
#include "./DipoleFromTable.h"
#include "./DarkNewsCrossSection.h"
#include "./DarkNewsDecay.h"
#include "./DISFromSpline.h"
#include "./DipoleDISFromSpline.h"
#include "./HNLFromSpline.h"
#include "./Decay.h"
#include "./NeutrissimoDecay.h"
#include "./InteractionCollection.h"
#include "./DummyCrossSection.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(interactions,m) {
    using namespace LI::interactions;

    register_CrossSection(m);
    register_Decay(m);
    register_DipoleFromTable(m);
    register_DarkNewsCrossSection(m);
    register_DarkNewsDecay(m);
    register_DISFromSpline(m);
    register_DipoleDISFromSpline(m);
    register_HNLFromSpline(m);
    register_NeutrissimoDecay(m);
    register_InteractionCollection(m);
    register_DummyCrossSection(m);
}
