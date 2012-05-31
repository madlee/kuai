#include <kuai/typedef.h>

#ifndef _KUAI_MOL_CHIRAL_CENTER_H_
#define _KUAI_MOL_CHIRAL_CENTER_H_

namespace kuai {

	class ChiralCenter {
	public:
		enum StereoType {
		   STEREO_TYPE_NONE        = 0,
		   STEREO_TYPE_DOUBLEBOND  = 1,
		   STEREO_TYPE_TETRAHEDRAL = 2,
		   STEREO_TYPE_ALLENE      = 3
		};
		enum StereoParity {
		   PARITY_NONE      = 0,
		   PARITY_ODD       = 1,  /* 'O' */
		   PARITY_EVEN      = 2,  /* 'E' */
		   PARITY_UNKNOWN   = 3,  /* 'U' */ /* (SEE ALSO READINCH.C)
												  USED IN: EXTRACT0DPARITIES, INCHITO_ATOM  */
		   PARITY_UNDEFINED = 4   /* '?' -- SHOULD NOT BE USED; HOWEVER, SEE NOTE ABOVE */
		};

	public:
		Index neighbor[4];		/* 4 atoms always */
		Index center_atom;		/* central tetrahedral atom or a central */
								/* atom of allene; otherwise INVALID_INDEX */
		StereoType type;
		StereoParity parity;    /* inchi_StereoParity0D: may be a combination of two parities: */
								/* ParityOfConnected | (ParityOfDisconnected << 3), see Note above */
	};

	typedef ChiralCenter* ChiralCenterPtr;
	typedef Array<ChiralCenter> ChiralCenterArray;
	typedef SharedPtr<ChiralCenterArray> ChiralCenterArrayPtr;
}

#endif
