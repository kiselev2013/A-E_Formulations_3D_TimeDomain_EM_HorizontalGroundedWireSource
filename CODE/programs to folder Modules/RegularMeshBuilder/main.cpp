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


#include "Logging.h"
#include "Processing.h"

int main(int argc, char **argv)
{
	char buf[2048];
	if (open_log("RegularMeshBuilder.log") != 0)
	{
		printf("Cound not open RegularMeshBuilder.log\n");
		return 1;
	}
	
	int status = Processing::MainProcedure();
	sprintf(buf, "Status = %d\n", status);
	write_to_log(buf);

	write_to_log("All done\n");
	fclose(log_file);

	return status;
}