/* Restraint to calculate a Harmonic restraint for particles if they move out of outer cylinder or if they move into the inner cylinder */

#ifndef IMPMICOS_Z_AXIAL_RESTRAINT_PROTEINWISE_H  // Header guard
#define IMPMICOS_Z_AXIAL_RESTRAINT_PROTEINWISE_H

#include <IMP/micos/micos_config.h>  // To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/Restraint.h>


IMPMICOS_BEGIN_NAMESPACE

class IMPMICOSEXPORT ZAxialRestraintProteinwise : public IMP::Restraint {
    IMP::ParticlesTemp plist_; // All the particles to which to apply this to (individual scores are summed)

    double lower_bound_;
    double upper_bound_;
    double sigma_;  // The sigma for the harmonic

    public:
        ZAxialRestraintProteinwise(IMP::ParticlesTemp plist, double lower_bound, double upper_bound, double sigma);

        // unprotected_evaluate calculates the score
        // do_get_inputs returns the particles for which the score was calculated
        virtual double unprotected_evaluate(IMP::Particle* p, IMP::DerivativeAccumulator* accum) const override;
        //IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
        virtual IMP::ModelObjectsTemp do_get_inputs() const override;
        IMP_OBJECT_METHODS(ZAxialRestraintProteinwise);  //add the usual IMP object methods

};

IMPMICOS_END_NAMESPACE

#endif
