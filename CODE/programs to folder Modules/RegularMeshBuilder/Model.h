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

#include <algorithm>

namespace Model
{
	using namespace std;

	// Gate for determining if materials properties are equal or not
#define MaterialEqualEps 1e-7

// Names of material properties
	enum PhysicalValueName
	{
		SigmaN,
		SigmaH,
		SigmaV,
		Eps,
		Mu
	};

	// Material property class
	struct PhysicalValue
	{
	public:
		double constantValue;

		PhysicalValue::PhysicalValue()
		{
			constantValue = 0.0;
		}
		PhysicalValue::PhysicalValue(const PhysicalValue &value)
		{
			constantValue = value.constantValue;
		}

		// Check if two constant values are equal
		inline bool ConstantEqual(double c1, double c2)
		{
			if (fabs(c1) < 1e-30 && fabs(c2) < 1e-30)
				return true;
			if (fabs(1.0 - min(c2, c1) / max(c2, c1)) < MaterialEqualEps)
				return true;
			return false;
		}

		// Check if two property value and this are equal
		inline bool PhysicalValue::Equal(PhysicalValue &value, bool compareTensors = false)
		{
			return ConstantEqual(constantValue, value.constantValue);
		}

	};

	// Material class
	struct Material
	{
	public:

		// Material properties:
		PhysicalValue SigmaH;
		PhysicalValue SigmaV;
		PhysicalValue SigmaN;
		PhysicalValue Eps;
		PhysicalValue Mu;
		int numberInNormalTask = 0;

		Material::Material()
		{
			Mu.constantValue = 1.0;
			Eps.constantValue = 0.0;
		}
		Material::Material(const Material &material)
		{
			CopyPhysicals(material);
		}

		// Copy material properties from another material to this
		void Material::CopyPhysicals(const Material &material)
		{
			this->SigmaH = material.SigmaH;
			this->SigmaV = material.SigmaV;
			this->SigmaN = material.SigmaN;
			this->Eps = material.Eps;
			this->Mu = material.Mu;
			this->numberInNormalTask = material.numberInNormalTask;
		}

		// Check if another material is equal to this
		bool Material::Equal(Material &material)
		{
			if (!material.SigmaH.Equal(SigmaH)) return false;
			if (!material.SigmaV.Equal(SigmaV)) return false;
			if (!material.SigmaN.Equal(SigmaN)) return false;
			if (!material.Eps.Equal(Eps)) return false;
			if (!material.Mu.Equal(Mu)) return false;

			return true;
		}

		PhysicalValue& Material::GetValue(PhysicalValueName name)
		{
			switch (name)
			{
			case PhysicalValueName::SigmaN:
				return SigmaN;
			case PhysicalValueName::SigmaH:
				return SigmaH;
			case PhysicalValueName::SigmaV:
				return SigmaV;
			case PhysicalValueName::Eps:
				return Eps;
			case PhysicalValueName::Mu:
				return Mu;
			}
		}
	};

	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));

	}

	// Box class for different purposes
	template <typename T>
	struct GeoBox
	{
	public:

		T coordinates[6]; // x0, x1, y0, y1, z0, z1
		int material = 0; // Material number
		Material materialValue; // Matrial value

		GeoBox()
		{

		}
	};

	// Layer class
	template <typename T>
	struct GeoLayer
	{
	public:

		T top; // Upper layer coordinate
		T bottom; // Lower layer coordinate
		int material = 0;  // Material number
		Material materialValue; // Material value
	};

	// Receiver class
	struct Receiver
	{
	public:

		double coordinates[3]; // x, y, z
	};

	struct Generator
	{
	public:

		double coordinates[3]; // x, y, z
	};

	// Create material for air
	void CreateAirMaterial(Material &airMaterial);
}