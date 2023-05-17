import cPickle, base64
try:
	from SimpleSession.versions.v65 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 16, 42360])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v65 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwJOfYdVCWJhbGxTY2FsZXEDSwJHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLAkc/8AAAAAAAAH2HVQVjb2xvcnEFSwJLAH2HVQpyaWJib25UeXBlcQZLAksAfYdVCnN0aWNrU2NhbGVxB0sCRz/wAAAAAAAAfYdVDG1tQ0lGSGVhZGVyc3EIXXEJKE5OZVUMYXJvbWF0aWNNb2RlcQpLAksBfYdVCnZkd0RlbnNpdHlxC0sCR0AUAAAAAAAAfYdVBmhpZGRlbnEMSwKJfYdVDWFyb21hdGljQ29sb3JxDUsCTn2HVQ9yaWJib25TbW9vdGhpbmdxDksCSwB9h1UJYXV0b2NoYWlucQ9LAoh9h1UKcGRiVmVyc2lvbnEQSwJLAH2HVQhvcHRpb25hbHERfXESVQhvcGVuZWRBc3ETiIlLAihVOS9ob21lL211c2thYW4vRG9jdW1lbnRzL21vZGVsc19yZXN1bHRzL3J1bl8xMC9ybWZzLzAucm1mM3EUTk5LAXRxFX2Hh3NVD2xvd2VyQ2FzZUNoYWluc3EWSwKJfYdVCWxpbmVXaWR0aHEXSwJHP/AAAAAAAAB9h1UPcmVzaWR1ZUxhYmVsUG9zcRhLAksAfYdVBG5hbWVxGUsCWEsAAAAvaG9tZS9tdXNrYWFuL0RvY3VtZW50cy9tb2RlbHNfcmVzdWx0cy9ydW5fMTAvcm1mcy8wLnJtZjMgLSBib3VuZHMgLSAwLnJtZjN9cRpYPgAAAC9ob21lL211c2thYW4vRG9jdW1lbnRzL21vZGVsc19yZXN1bHRzL3J1bl8xMC9ybWZzLzAucm1mMyAtIDEzXXEbSwBhc4dVD2Fyb21hdGljRGlzcGxheXEcSwKJfYdVD3JpYmJvblN0aWZmbmVzc3EdSwJHP+mZmZmZmZp9h1UKcGRiSGVhZGVyc3EeXXEfKH1xIH1xIWVVA2lkc3EiSwJLAEsBhn1xI0sASwKGXXEkSwFhc4dVDnN1cmZhY2VPcGFjaXR5cSVLAke/8AAAAAAAAH2HVRBhcm9tYXRpY0xpbmVUeXBlcSZLAksCfYdVFHJpYmJvbkhpZGVzTWFpbmNoYWlucSdLAoh9h1UHZGlzcGxheXEoSwKIfYd1Lg=='))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAktlVQEgfYdVC2ZpbGxEaXNwbGF5cQNLZYl9h1UEbmFtZXEES2VYAwAAAEJTUH1xBShYAwAAAENZU11xBihLAEsqZVgDAAAAR0xOXXEHSythWAMAAABBU1BdcQhLA2FYAwAAAFNFUl1xCShLE0sfSyhlWAMAAABWQUxdcQooSwVLBksRZVgDAAAATFlTXXELSwdhWAMAAABQUk9dcQxLGmFYAwAAAFRIUl1xDShLCksVZVgDAAAAUEhFXXEOKEsMSxJLFksXSx1lWAMAAABBTEFdcQ8oSwJLBEscSyZlWAMAAABISVNdcRBLLGFYAwAAAEdMWV1xEShLCUsLSw1LD0seSyBLIkskZVgDAAAASUxFXXESKEsISxBlWAMAAABMRVVdcRMoSwFLDksUSxtLI2VYAwAAAFRSUF1xFEsZYVgDAAAAQVNOXXEVSylhWAMAAABUWVJdcRZLJ2FYAwAAAE1FVF1xFyhLGEshSyVldYdVBWNoYWlucRhLZVgBAAAAIH1xGVgBAAAAQU5dcRpLAEsthnEbYYZzh1UOcmliYm9uRHJhd01vZGVxHEtlSwJ9h1UCc3NxHUtliYmGfYdVCG1vbGVjdWxlcR5LZUsBfXEfSwBOXXEgSwBLLYZxIWGGc4dVC3JpYmJvbkNvbG9ycSJLZUsBfXEjKEsCTl1xJChLLUsFhnElS0lLAYZxJktLSwOGcSdLY0sChnEoZYZLA05dcSkoSzJLF4ZxKktOSxWGcStlhksETl1xLEtKSwGGcS1hhnWHVQVsYWJlbHEuS2VYAAAAAH2HVQpsYWJlbENvbG9ycS9LZU59cTAoSwNOXXExKEsySxeGcTJLTksVhnEzZYZLBE5dcTRLSksBhnE1YYZ1h1UIZmlsbE1vZGVxNktlSwF9h1UFaXNIZXRxN0tliX2HVQtsYWJlbE9mZnNldHE4S2VOfYdVCHBvc2l0aW9ucTldcTooSw1LGIZxO0soSxWGcTxLAEs4hnE9ZVUNcmliYm9uRGlzcGxheXE+S2WJfYdVCG9wdGlvbmFscT99VQRzc0lkcUBLZUr/////fYd1Lg=='))
	atomInfo = cPickle.loads(base64.b64decode('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YWxxqX1xqihVDHNlcmlhbE51bWJlcnGriIhdcawoSv////9LAYZxrUr/////SwGGca5K/////0sBhnGvSv////9LAYZxsEr/////SwGGcbFK/////0sBhnGySv////9LAYZxs0r/////SwGGcbRK/////0sBhnG1Sv////9LAYZxtkr/////SwGGcbdK/////0sBhnG4Sv////9LAYZxuUr/////SwGGcbpK/////0sBhnG7Sv////9LAYZxvEr/////SwGGcb1K/////0sBhnG+Sv////9LAYZxv0r/////SwGGccBK/////0sBhnHBSv////9LAYZxwkr/////SwGGccNK/////0sBhnHESv////9LAYZxxUr/////SwGGccZK/////0sBhnHHSv////9LAYZxyEr/////SwGGcclK/////0sBhnHKSv////9LAYZxy0r/////SwGGccxK/////0sBhnHNSv////9LAYZxzkr/////SwGGcc9K/////0sBhnHQSv////9LAYZx0Ur/////SwGGcdJK/////0sBhnHTSv////9LAYZx1Er/////SwGGcdVK/////0sBhnHWSv////9LAYZx10r/////SwGGcdhK/////0sBhnHZSv////9LAYZx2kr/////SwGGcdtK/////0sBhnHcSv////9LAYZx3Ur/////SwGGcd5K/////0sBhnHfSv////9LAYZx4Er/////SwGGceFK/////0sBhnHiSv////9LAYZx40r/////SwGGceRlh1UHYmZhY3RvcnHliIlLOEcAAAAAAAAAAH2Hh1UJb2NjdXBhbmN5ceaIiUs4Rz/wAAAAAAAAfYeHdVUHZGlzcGxheXHnSziIfXHoiU5dcekoSwJLA4Zx6kseSwOGcetlhnOHdS4='))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECSwBOfYdVBWF0b21zcQNdVQVsYWJlbHEESwBOfYdVCGhhbGZib25kcQVLAE59h1UGcmFkaXVzcQZLAE59h1ULbGFiZWxPZmZzZXRxB0sATn2HVQhkcmF3TW9kZXEISwBOfYdVCG9wdGlvbmFscQl9VQdkaXNwbGF5cQpLAE59h3Uu'))
	crdInfo = cPickle.loads(base64.b64decode('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'))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (1, 0, {}), 'display': (1, True, {}), 'id': (1, 1, {}), 'vrmlString': ['.color white\n.cylinder 0 0 -200 0 0 200 12 open\n.color red\n.cylinder 0 0 -200 0 0 200 19 open\n'], 'name': (1, u'/home/muskaan/Documents/micos/modeling/two_cylinder_19,12.bild', {})}
	colors = {u'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), u'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), u'Rf': ((0.8, 0, 0.34902), 1, u'default'), u'Ra': ((0, 0.490196, 0), 1, u'default'), u'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), u'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), u'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), u'Be': ((0.760784, 1, 0), 1, u'default'), u'Ba': ((0, 0.788235, 0), 1, u'default'), u'Bh': ((0.878431, 0, 0.219608), 1, u'default'), u'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), u'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), u'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), u'H': ((1, 1, 1), 1, u'default'), u'P': ((1, 0.501961, 0), 1, u'default'), u'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), u'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), u'Gd': ((0.270588, 1, 0.780392), 1, u'default'), u'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'), u'Pr': ((0.85098, 1, 0.780392), 1, u'default'), u'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'), u'Pu': ((0, 0.419608, 1), 1, u'default'),
u'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), u'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), u'Pa': ((0, 0.631373, 1), 1, u'default'), u'Pd': ((0, 0.411765, 0.521569), 1, u'default'), u'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), u'Po': ((0.670588, 0.360784, 0), 1, u'default'), u'Pm': ((0.639216, 1, 0.780392), 1, u'default'), u'Hs': ((0.901961, 0, 0.180392), 1, u'default'), u'Ho': ((0, 1, 0.611765), 1, u'default'), u'Hf': ((0.301961, 0.760784, 1), 1, u'default'), u'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), u'He': ((0.85098, 1, 1), 1, u'default'), u'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), u'Mg': ((0.541176, 1, 0), 1, u'default'), u'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), u'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'), u'O': ((1, 0.0509804, 0.0509804), 1, u'default'), u'Mt': ((0.921569, 0, 0.14902), 1, u'default'), u'S': ((1, 1, 0.188235), 1, u'default'), u'W': ((0.129412, 0.580392, 0.839216), 1, u'default'), u'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'), u'Eu': ((0.380392, 1, 0.780392), 1, u'default'),
u'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), u'Er': ((0, 0.901961, 0.458824), 1, u'default'), u'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), u'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), u'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), u'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), u'Nd': ((0.780392, 1, 0.780392), 1, u'default'), u'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), u'Np': ((0, 0.501961, 1), 1, u'default'), u'Fr': ((0.258824, 0, 0.4), 1, u'default'), u'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), u'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), u'B': ((1, 0.709804, 0.709804), 1, u'default'), u'F': ((0.564706, 0.878431, 0.313725), 1, u'default'), u'Sr': ((0, 1, 0), 1, u'default'), u'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), u'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), u'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'), u'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'), u'Sm': ((0.560784, 1, 0.780392), 1, u'default'), u'V': ((0.65098, 0.65098, 0.670588), 1, u'default'),
u'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), u'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), u'Sg': ((0.85098, 0, 0.270588), 1, u'default'), u'Se': ((1, 0.631373, 0), 1, u'default'), u'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), u'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), u'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), u'Ca': ((0.239216, 1, 0), 1, u'default'), u'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), u'Ce': ((1, 1, 0.780392), 1, u'default'), u'Cd': ((1, 0.85098, 0.560784), 1, u'default'), u'Lu': ((0, 0.670588, 0.141176), 1, u'default'), u'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'), u'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), u'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'), u'La': ((0.439216, 0.831373, 1), 1, u'default'), u'Li': ((0.8, 0.501961, 1), 1, u'default'), u'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'), u'Tm': ((0, 0.831373, 0.321569), 1, u'default'), u'Lr': ((0.780392, 0, 0.4), 1, u'default'), u'Th': ((0, 0.729412, 1), 1, u'default'), u'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'),
u'Te': ((0.831373, 0.478431, 0), 1, u'default'), u'Tb': ((0.188235, 1, 0.780392), 1, u'default'), u'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), u'Ta': ((0.301961, 0.65098, 1), 1, u'default'), u'Yb': ((0, 0.74902, 0.219608), 1, u'default'), u'Db': ((0.819608, 0, 0.309804), 1, u'default'), u'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'), u'Dy': ((0.121569, 1, 0.780392), 1, u'default'), u'I': ((0.580392, 0, 0.580392), 1, u'default'), u'U': ((0, 0.560784, 1), 1, u'default'), u'Y': ((0.580392, 1, 1), 1, u'default'), u'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), u'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), u'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'), u'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), u'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), u'As': ((0.741176, 0.501961, 0.890196), 1, u'default'), u'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'), u'Au': ((1, 0.819608, 0.137255), 1, u'default'), u'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), u'In': ((0.65098, 0.458824, 0.45098), 1, u'default'),
u'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default')}
	materials = {u'default': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (0, None, {}), 'atoms': [], 'label': (0, None, {}), 'halfbond': (0, None, {}), 'labelColor': (0, None, {}), 'labelOffset': (0, None, {}), 'drawMode': (0, None, {}), 'display': (0, None, {})}], 'lineType': (1, 2, {}), 'color': (1, 4, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (5, (u'', (1, 0, 0, 1)), {(u'green', (0, 1, 0, 1)): [3], (u'gray', (0.745, 0.745, 0.745, 1)): [0], (u'yellow', (1, 1, 0, 1)): [4]})
	viewerInfo = {'cameraAttrs': {'center': (2.5541044473648, 11.693176388741, 26.110602855682), 'fieldOfView': 29.495693998264, 'nearFar': (204.07896762017, -232.96341730161), 'ortho': False, 'eyeSeparation': 50.8, 'focal': 26.110602855682}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 55.523655742771, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': False, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [(None, [(0.94117647058824, 0.94117647058824, 0.94117647058824, 1), (0.74117647058824, 0.74117647058824, 0.74117647058824, 1), (0.38823529411765, 0.38823529411765, 0.38823529411765, 1)], 1), 1, 0, 0], 'depthCue': False, 'highlight': 0, 'scaleFactor': 0.50520648544392, 'angleDependentTransparency': False, 'backgroundMethod': 0}, 'viewerHL': 3, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': None}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v65 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [ ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_cap_attributes():
 cap_attributes = \
  {
   'cap_attributes': [ ],
   'cap_color': None,
   'cap_offset': 0.01,
   'class': 'Caps_State',
   'default_cap_offset': 0.01,
   'mesh_style': False,
   'shown': True,
   'subdivision_factor': 1.0,
   'version': 1,
  }
 import SurfaceCap.session
 SurfaceCap.session.restore_cap_attributes(cap_attributes)
registerAfterModelsCB(restore_cap_attributes)

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = []
	userXSections = []
	userResidueClasses = []
	residueData = [(2, 'Chimera default', 'rounded', u'unknown'), (3, 'Chimera default', 'rounded', u'unknown'), (4, 'Chimera default', 'rounded', u'unknown'), (5, 'Chimera default', 'rounded', u'unknown'), (6, 'Chimera default', 'rounded', u'unknown'), (7, 'Chimera default', 'rounded', u'unknown'), (8, 'Chimera default', 'rounded', u'unknown'), (9, 'Chimera default', 'rounded', u'unknown'), (10, 'Chimera default', 'rounded', u'unknown'), (11, 'Chimera default', 'rounded', u'unknown'), (12, 'Chimera default', 'rounded', u'unknown'), (13, 'Chimera default', 'rounded', u'unknown'), (14, 'Chimera default', 'rounded', u'unknown'), (15, 'Chimera default', 'rounded', u'unknown'), (16, 'Chimera default', 'rounded', u'unknown'), (17, 'Chimera default', 'rounded', u'unknown'), (18, 'Chimera default', 'rounded', u'unknown'), (19, 'Chimera default', 'rounded', u'unknown'), (20, 'Chimera default', 'rounded', u'unknown'), (21, 'Chimera default', 'rounded', u'unknown'), (22, 'Chimera default', 'rounded', u'unknown'), (23, 'Chimera default', 'rounded', u'unknown'), (24, 'Chimera default', 'rounded', u'unknown'),
(25, 'Chimera default', 'rounded', u'unknown'), (26, 'Chimera default', 'rounded', u'unknown'), (27, 'Chimera default', 'rounded', u'unknown'), (28, 'Chimera default', 'rounded', u'unknown'), (29, 'Chimera default', 'rounded', u'unknown'), (30, 'Chimera default', 'rounded', u'unknown'), (31, 'Chimera default', 'rounded', u'unknown'), (32, 'Chimera default', 'rounded', u'unknown'), (33, 'Chimera default', 'rounded', u'unknown'), (34, 'Chimera default', 'rounded', u'unknown'), (35, 'Chimera default', 'rounded', u'unknown'), (36, 'Chimera default', 'rounded', u'unknown'), (37, 'Chimera default', 'rounded', u'unknown'), (38, 'Chimera default', 'rounded', u'unknown'), (39, 'Chimera default', 'rounded', u'unknown'), (40, 'Chimera default', 'rounded', u'unknown'), (41, 'Chimera default', 'rounded', u'unknown'), (42, 'Chimera default', 'rounded', u'unknown'), (43, 'Chimera default', 'rounded', u'unknown'), (44, 'Chimera default', 'rounded', u'unknown'), (45, 'Chimera default', 'rounded', u'unknown'), (46, 'Chimera default', 'rounded', u'unknown'), (47, 'Chimera default', 'rounded', u'unknown'),
(48, 'Chimera default', 'rounded', u'unknown'), (49, 'Chimera default', 'rounded', u'unknown'), (50, 'Chimera default', 'rounded', u'unknown'), (51, 'Chimera default', 'rounded', u'unknown'), (52, 'Chimera default', 'rounded', u'unknown'), (53, 'Chimera default', 'rounded', u'unknown'), (54, 'Chimera default', 'rounded', u'unknown'), (55, 'Chimera default', 'rounded', u'unknown'), (56, 'Chimera default', 'rounded', u'unknown'), (57, 'Chimera default', 'rounded', u'unknown'), (58, 'Chimera default', 'rounded', u'unknown'), (59, 'Chimera default', 'rounded', u'unknown'), (60, 'Chimera default', 'rounded', u'unknown'), (61, 'Chimera default', 'rounded', u'unknown'), (62, 'Chimera default', 'rounded', u'unknown'), (63, 'Chimera default', 'rounded', u'unknown'), (64, 'Chimera default', 'rounded', u'unknown'), (65, 'Chimera default', 'rounded', u'unknown'), (66, 'Chimera default', 'rounded', u'unknown'), (67, 'Chimera default', 'rounded', u'unknown'), (68, 'Chimera default', 'rounded', u'unknown'), (69, 'Chimera default', 'rounded', u'unknown'), (70, 'Chimera default', 'rounded', u'unknown'),
(71, 'Chimera default', 'rounded', u'unknown'), (72, 'Chimera default', 'rounded', u'unknown'), (73, 'Chimera default', 'rounded', u'unknown'), (74, 'Chimera default', 'rounded', u'unknown'), (75, 'Chimera default', 'rounded', u'unknown'), (76, 'Chimera default', 'rounded', u'unknown'), (77, 'Chimera default', 'rounded', u'unknown'), (78, 'Chimera default', 'rounded', u'unknown'), (79, 'Chimera default', 'rounded', u'unknown'), (80, 'Chimera default', 'rounded', u'unknown'), (81, 'Chimera default', 'rounded', u'unknown'), (82, 'Chimera default', 'rounded', u'unknown'), (83, 'Chimera default', 'rounded', u'unknown'), (84, 'Chimera default', 'rounded', u'unknown'), (85, 'Chimera default', 'rounded', u'unknown'), (86, 'Chimera default', 'rounded', u'unknown'), (87, 'Chimera default', 'rounded', u'unknown'), (88, 'Chimera default', 'rounded', u'unknown'), (89, 'Chimera default', 'rounded', u'unknown'), (90, 'Chimera default', 'rounded', u'unknown'), (91, 'Chimera default', 'rounded', u'unknown'), (92, 'Chimera default', 'rounded', u'unknown'), (93, 'Chimera default', 'rounded', u'unknown'),
(94, 'Chimera default', 'rounded', u'unknown'), (95, 'Chimera default', 'rounded', u'unknown'), (96, 'Chimera default', 'rounded', u'unknown'), (97, 'Chimera default', 'rounded', u'unknown'), (98, 'Chimera default', 'rounded', u'unknown'), (99, 'Chimera default', 'rounded', u'unknown'), (100, 'Chimera default', 'rounded', u'unknown'), (101, 'Chimera default', 'rounded', u'unknown'), (102, 'Chimera default', 'rounded', u'unknown')]
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4='
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.35740674433659325, 0.6604015517481454, -0.6604015517481455), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.35740674433659325, 0.6604015517481454, 0.6604015517481455), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.25056280708573153, 0.25056280708573153, 0.9351131265310293), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
	from SimpleSession.versions.v65 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  [132]
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 4 }
	windowSize = (1405, 703)
	xformMap = {0: (((-0.10179108086136, 0.10335204600509, 0.98942252371959), 6.8113204903344), (-15.318066778727, -1.9307538841621, -14.442224840719), True), 1: (((-0.10179108086136, 0.10335204600509, 0.98942252371959), 6.8113204903344), (-15.318066778727, -1.9307538841621, -14.442224840719), True), 159: (((-0.10179108086136, 0.10335204600509, 0.98942252371959), 6.8113204903344), (-15.318066778727, -1.9307538841621, -14.442224840719), True)}
	fontInfo = {'face': ('Sans Serif', 'Normal', 16)}
	clipPlaneInfo = {}
	silhouettes = {0: True, 1: True, 160: True, 159: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

