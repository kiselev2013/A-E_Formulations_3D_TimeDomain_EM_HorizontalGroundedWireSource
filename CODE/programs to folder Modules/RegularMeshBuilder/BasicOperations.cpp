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

#include "BasicOperations.h"
#include "Logging.h"

namespace BasicOperations
{
	using namespace std;
	// Find index of interval containing value in v
	int BinarySearchInVectorSorted(vector<double> &v, double value)
	{
		int b[2] = { 0, v.size() - 1 };
		int m;

		if (value < v.front()) return 0;
		if (value > v.back()) return v.size() - 1;

		while (b[1] - b[0] > 1)
		{
			m = (b[1] + b[0]) / 2;
			if (value > v[m])
				b[0] = m;
			else
				b[1] = m;
		}

		if (b[0] == b[1])
			return b[0];

		if (fabs(v[b[0]] - value) < fabs(v[b[1]] - value))
			return b[0];

		return b[1];
	}

	// Sort coordinates
	int InsertCoordinateIntoSortedArray(vector<double> &arr, double value, double eps)
	{
		try
		{
			if (value <= arr.front())
				return 0;
			if (value >= arr.back())
				return 0;

			int	foundIndex = BinarySearchInVectorSorted(arr, value);
			int neighbourIndex = arr[foundIndex] > value && foundIndex != 0 ? foundIndex - 1 : foundIndex + 1;
			double step = fabs(arr[foundIndex] - arr[neighbourIndex]);
			if (fabs(arr[foundIndex] - value) / step > eps)
			{
				arr.insert(arr.begin() + (arr[foundIndex] > value ? foundIndex : foundIndex + 1), value);
			}
			else
				arr[foundIndex] = value;

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : InsertCoordinateIntoSortedArray : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Round double value with decrement
	double RoundDouble_(double value, double decrement)
	{
		int n = fabs(value / decrement);
		return value < 0 ? -(n + 1) * decrement : n * decrement;
	}
}