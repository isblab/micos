#include <cmath>
#include <algorithm>
#include <IMP/Model.h>
#include <IMP/Particle.h>
#include <IMP/SingletonModifier.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/XYZR.h>
#include <IMP/internal/container_helpers.h>
#include <IMP/internal/StaticListContainer.h>
#include <IMP/micos/MembraneExclusionRestraint.h>


IMPMICOS_BEGIN_NAMESPACE

namespace {
IMP::internal::StaticListContainer<SingletonContainer> *get_list(
    SingletonContainer *sc) {
  IMP::internal::StaticListContainer<SingletonContainer> *ret =
      dynamic_cast<
          IMP::internal::StaticListContainer<SingletonContainer> *>(
          sc);
  if (!ret) {
    IMP_THROW("Can only use the set and add methods when no container"
                  << " was passed on construction of ConnectivityRestraint.",
              ValueException);
  }
  return ret;
}
}

/*##################RESTRAINT SETUP ###########################*/


MembraneExclusionRestraint::MembraneExclusionRestraint(Model *m,
    SingletonContainerAdaptor sc,
    double R, double r, double thickness, double sigma)
  : Restraint(m, "MembraneExclusionRestraint %1%")
  , R_(R)
  , r_(r)
  , sigma_(sigma)
  , thickness_(thickness/2)
{
  sc_ = sc;
}

MembraneExclusionRestraint::MembraneExclusionRestraint(Model *m,
    double R, double r, double thickness, double sigma)
  : Restraint(m, "MembraneExclusionRestraint %1%")
  , R_(R)
  , r_(r)
  , sigma_(sigma)
  , thickness_(thickness/2)
{
}

void MembraneExclusionRestraint::set_particles(const ParticlesTemp &ps) {
  if (!sc_ && !ps.empty()) {
    sc_ = new IMP::internal::StaticListContainer<SingletonContainer>(
        ps[0]->get_model(), "Membrane exclusion list");
  }
  get_list(sc_)->set(IMP::internal::get_index(ps));
}

void MembraneExclusionRestraint::add_particles(const ParticlesTemp &ps) {
  if (!sc_ && !ps.empty()) {
    sc_ = new IMP::internal::StaticListContainer<SingletonContainer>(
        ps[0]->get_model(), "Membrane exclusion list");
  }
  get_list(sc_)->add(IMP::internal::get_index(ps));
}

void MembraneExclusionRestraint::add_particle(Particle *ps) {
  if (!sc_) {
    sc_ = new IMP::internal::StaticListContainer<SingletonContainer>(
        ps->get_model(), "Membrane exclusion list");
  }
  get_list(sc_)->add(IMP::internal::get_index(ps));
}

double
MembraneExclusionRestraint::unprotected_evaluate(DerivativeAccumulator *) const
{
  IMP_CHECK_OBJECT(sc_.get());
  double v = 0.0, x, y, z;
  core::XYZ i_current;
  IMP::ParticlesTemp all_particles = sc_->get();
  double excluded_radius = r_ + thickness_;
  double membrane_boundary = (R_ - excluded_radius) * (R_ - excluded_radius);

  for (unsigned int i = 0; i < all_particles.size(); ++i )
  {
    i_current = core::XYZ(all_particles[i]);
    x = i_current.get_coordinate(0);
    y = i_current.get_coordinate(1);
    z = i_current.get_coordinate(2);

    if ( std::fabs(z) > excluded_radius )
      continue;

    if ( (x*x + y*y) < membrane_boundary )
      continue;

    std::pair<double, algebra::Vector3D> dist = half_torus_distance(x, y, z, R_, r_);
    if ( dist.first < thickness_ )
    {
      v += (thickness_ - dist.first) * (thickness_ - dist.first);
      /*
      if ( accum )
      {
        all_particles[i]->get_model()->add_to_coordinate_derivatives(IMP::internal::get_index(all_particles[i]),
            dist.second*2*dist.first/sigma_, *accum);
      }
      */
    }
  }
  return v/sigma_;
}

ModelObjectsTemp MembraneExclusionRestraint::do_get_inputs() const {
  if ( !sc_ )
    return ModelObjectsTemp();
  ParticleIndexes all = sc_->get_all_possible_indexes();
  return IMP::get_particles(get_model(), all);
}


IMPMICOS_END_NAMESPACE
