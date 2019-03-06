/*
 * Options.h
 *
 *  Created on: Jun 21, 2017
 *      Author: cjustin
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>
#include <stdint.h>

typedef uint32_t BarcodeID;
typedef uint64_t Minimizer;
typedef uint32_t Count;

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
	extern unsigned minN;
	extern unsigned threads;
}

#endif /* OPTIONS_H_ */
