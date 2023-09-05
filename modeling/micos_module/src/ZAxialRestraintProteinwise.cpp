#include <cmath>
#include <IMP/core/XYZ.h>
#include <IMP/micos/ZAxialRestraintProteinwise.h>

IMPMICOS_BEGIN_NAMESPACE

// tries to keep the particles above CJ
ZAxialRestraintProteinwise::ZAxialRestraintProteinwise(IMP::ParticlesTemp plist, double upper_bound, double lower_bound, double sigma) :

		Restraint(plist[0]->get_model(), "ZAxialRestraintProteinwise %1%"),
		plist_(plist),
		upper_bound_(upper_bound),
		lower_bound_(lower_bound),
		sigma_(sigma){}


double ZAxialRestraintProteinwise::unprotected_evaluate(IMP::Particle* p, IMP::DerivativeAccumulator* accum) const {
	double score = 0;
	for (unsigned int i=0; i < plist_.size(); i++){
		score += IMP::core::XYZ(plist_[i]).get_coordinate(2);
		        }

	double avg = score/plist_.size()

	if (avg<lower_bound_) {
		double deviation = avg*avg + lower_bound_*lower_bound_ - 2*avg*lower_bound_;
		return fabs(deviation)/sigma_;
			}

	if (avg>upper_bound_) {
		double deviation = avg*avg + upper_bound_*upper_bound_ - 2*avg*upper_bound_;
		return fabs(deviation)/sigma_;
		}

	else{
		return 0;
	}
}

IMP::ModelObjectsTemp ZAxialRestraintProteinwise::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
