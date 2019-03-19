/*
 * Options.h
 *
 *  Created on: Jun 21, 2017
 *      Author: cjustin
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <stdint.h>

using BarcodeID = uint32_t;
using Minimizer = uint64_t;
using Count = uint16_t;

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
	extern unsigned minN;
	extern unsigned threads;
}

#endif /* OPTIONS_H_ */
