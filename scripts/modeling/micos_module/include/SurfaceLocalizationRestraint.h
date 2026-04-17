/* Restraint to calculate a Harmonic restraint for particles if they are away from the surface of inner cylinder*/

#ifndef IMPMICOS_SURFACE_LOCALIZATION_RESTRAINT_H  // Header guard
#define IMPMICOS_SURFACE_LOCALIZATION_RESTRAINT_H

#include <IMP/micos/micos_config.h>  // To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/Restraint.h>
#include "IMP/Model.h"


IMPMICOS_BEGIN_NAMESPACE

class IMPMICOSEXPORT SurfaceLocalizationRestraint : public IMP::Restraint {
    IMP::ParticlesTemp plist_; // All the particles to which to apply this to (individual scores are summed)
    double rsq_;  // inner radius squared
    double sigma_;  // The sigma for the harmonic
    double min_dist_sq;

    public:
        SurfaceLocalizationRestraint(IMP::ParticlesTemp plist, double r, double sigma, double allowed_dist);

        // unprotected_evaluate calculates the score
        // do_get_inputs returns the particles for which the score was calculated
        virtual double unprotected_evaluate(IMP::DerivativeAccumulator* accum) const override;
        //IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
        virtual IMP::ModelObjectsTemp do_get_inputs() const override;
        IMP_OBJECT_METHODS(SurfaceLocalizationRestraint);  //add the usual IMP object methods

    private:
 	double getDistance(IMP::Particle* p) const;
};

IMPMICOS_END_NAMESPACE

#endif
