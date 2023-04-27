#include <cmath>
#include <algorithm>
#include <IMP/Model.h>
#include <IMP/Particle.h>
#include <IMP/SingletonModifier.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/XYZR.h>
#include <IMP/internal/container_helpers.h>
#include <IMP/internal/StaticListContainer.h>
#include <IMP/micos/Restraint.h>

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


/*#####################################################
# Restraints setup - MembraneSurfaceLocationRestraint
#####################################################*/
MembraneSurfaceLocationRestraint::MembraneSurfaceLocationRestraint(Model *m,
    SingletonContainerAdaptor sc,
    double R, double r, double thickness, double sigma)
  : Restraint(m, "MembraneSurfaceLocationRestraint %1%")
  , R_(R)
  , r_(r)
  , sigma_(sigma)
  , thickness_(thickness/2)
{
  sc_ = sc;
}

MembraneSurfaceLocationRestraint::MembraneSurfaceLocationRestraint(Model *m,
    double R, double r, double thickness, double sigma)
  : Restraint(m, "MembraneSurfaceLocationRestraint %1%")
  , R_(R)
  , r_(r)
  , sigma_(sigma)
  , thickness_(thickness/2)
{
}

void MembraneSurfaceLocationRestraint::set_particles(const ParticlesTemp &ps) {
  if (!sc_ && !ps.empty()) {
    sc_ = new IMP::internal::StaticListContainer<SingletonContainer>(
        ps[0]->get_model(), "Surface location list");
  }
  get_list(sc_)->set(IMP::internal::get_index(ps));
}

void MembraneSurfaceLocationRestraint::add_particles(const ParticlesTemp &ps) {
  if (!sc_ && !ps.empty()) {
    sc_ = new IMP::internal::StaticListContainer<SingletonContainer>(
        ps[0]->get_model(), "Surface location list");
  }
  get_list(sc_)->add(IMP::internal::get_index(ps));
}

void MembraneSurfaceLocationRestraint::add_particle(Particle *ps) {
  if (!sc_) {
    sc_ = new IMP::internal::StaticListContainer<SingletonContainer>(
        ps->get_model(), "Surface location list");
  }
  get_list(sc_)->add(IMP::internal::get_index(ps));
}

double
MembraneSurfaceLocationRestraint::unprotected_evaluate(DerivativeAccumulator *) const
{
  IMP_CHECK_OBJECT(sc_.get());
  double v = 0;
  IMP::ParticlesTemp all_particles = sc_->get();
  for (unsigned int i = 0; i < all_particles.size(); ++i )
  {
    core::XYZ i_current = core::XYZ(all_particles[i]);
    /*
    double z = i_current.get_coordinate(2);
    if ( z > r_ )
      z = r_;
    else if ( z < -r_ )
      z = -r_;
    std::pair<double, algebra::Vector3D> dist = half_torus_distance(i_current.get_coordinate(0),
        i_current.get_coordinate(1), z, R_, r_);
    */
    std::pair<double, algebra::Vector3D> dist = half_torus_distance(i_current.get_coordinate(0),
        i_current.get_coordinate(1), i_current.get_coordinate(2), R_, r_);
    if ( std::fabs(dist.first) > thickness_ )
    {
      v += dist.first*dist.first;
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

ModelObjectsTemp MembraneSurfaceLocationRestraint::do_get_inputs() const {
  if ( !sc_ )
    return ModelObjectsTemp();
  ParticleIndexes all = sc_->get_all_possible_indexes();
  return IMP::get_particles(get_model(), all);
}



IMPMICOS_END_NAMESPACE
