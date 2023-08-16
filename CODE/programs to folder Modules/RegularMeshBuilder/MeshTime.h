/**
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *
 * The file contains basic utility functions and simple operations
 *
 *  Written by Ph.D. Dmitry S. Kiselev
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  harlequin_00@mail.ru
  *  Version 2.0 16 January, 2023
*/

#pragma once


namespace MeshTime
{
	struct Settings
	{
	public:
		double timeStep = 1e-6;    // Time mesh initial steps
		double timeSparse = 1.1;   // Time mesh sparce coefficient
		double timeLast = 1.43e+4; // last time in time mesh
	};


}

