/*
 * GENERAL REMARKS
 *  
 *  This code is freely available under the following conditions:
 *  
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *  
 *  			TDEMLineCalc
 *  This file contains code for read-write recevires data functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once
#define nul 1e-10

// The class contains the node data of a 3D grid
class Point3D
{
private:
protected:
	double xcoord, ycoord, zcoord; //!<  
public:
	Point3D()
	{ 
		xcoord=ycoord=zcoord=0.; 
	}
	Point3D(const double &xc, const double &yc, const double &zc) 
	{ 
		xcoord=xc; ycoord=yc; zcoord=zc; 
	}
	virtual ~Point3D() {}
	double& x() { return xcoord; }
	double& y() { return ycoord; }
	double& z() { return zcoord; }
	double getx() const { return xcoord; }
	double gety() const { return ycoord; }
	double getz() const { return zcoord; }
	Point3D operator=(const Point3D &p)
	{
		xcoord=p.getx();
		ycoord=p.gety();
		zcoord=p.getz();
		return *this;
	}
	Point3D operator=(const double &v)
	{
		xcoord=v;
		ycoord=v;
		zcoord=v;
		return *this;
	}
	Point3D operator+(const Point3D &p) const
	{ 
		return Point3D(xcoord+p.getx(), ycoord+p.gety(), zcoord+p.getz());
	}
	Point3D operator-(const Point3D &p) const
	{ 
		return Point3D(xcoord-p.getx(), ycoord-p.gety(), zcoord-p.getz());
	}
	Point3D operator+=(const Point3D &p)
	{ 
		xcoord+=p.getx(); ycoord+=p.gety(); zcoord+=p.getz();
		return *this;
	}
	Point3D operator-=(const Point3D &p)
	{ 
		xcoord-=p.getx(); ycoord-=p.gety(); zcoord-=p.getz();
		return *this;
	}
	Point3D operator*(const double &a) const
	{
		return Point3D(xcoord*a, ycoord*a, zcoord*a);
	}
	Point3D operator*=(const double &a)
	{
		xcoord*=a; ycoord*=a; zcoord*=a;
		return *this;
	}
	double& operator[](const int &i)
	{
		if (i==0) return xcoord;
		if (i==1) return ycoord;
		return zcoord;
	}
	const double& operator[](const int &i) const
	{
		if (i==0) return xcoord;
		if (i==1) return ycoord;
		return zcoord;
	}
	bool operator==(const Point3D &P) const
	{
		return (fabs(xcoord-P.getx())<nul&&
			fabs(ycoord-P.gety())<nul&&
			fabs(zcoord-P.getz())<nul);
	}
	friend ofstream& operator<<(ofstream &file, const Point3D &P) 
	{
		file << P.xcoord << " " << P.ycoord << " " << P.zcoord;
		return file;
	}
	friend ifstream& operator>>(ifstream &file, Point3D &P)
	{
		file >> P.xcoord >> P.ycoord >> P.zcoord;
		return file;
	}

	friend Point3D minp(const Point3D &P1, const Point3D &P2)
	{
		return Point3D(min(P1.xcoord, P2.xcoord),
			min(P1.ycoord, P2.ycoord),
			min(P1.zcoord, P2.zcoord));
	}
	friend Point3D maxp(const Point3D &P1, const Point3D &P2)
	{
		return Point3D(max(P1.xcoord, P2.xcoord),
			max(P1.ycoord, P2.ycoord),
			max(P1.zcoord, P2.zcoord));
	}

	void GetCoords(float *m)
	{
		m[0]=static_cast<float>(xcoord);
		m[1]=static_cast<float>(ycoord);
		m[2]=static_cast<float>(zcoord);
	}
};

class Vector: public Point3D
{
private:
protected:
public:
	Vector(): Point3D() {}
	Vector(const double &xc, const double &yc, const double &zc): Point3D(xc, yc, zc) {} 
	Vector(const Point3D &P) { xcoord=P.getx(); ycoord=P.gety(); zcoord=P.getz(); } 
	virtual ~Vector() {}
	Vector operator=(const Point3D &p)
	{
		dynamic_cast<Point3D*>(this)->operator=(p);
		return *this;
	}
	void operator/=(const double &d)
	{
		xcoord/=d; ycoord/=d; zcoord/=d; 
	}
	void operator*=(const double &d)
	{
		xcoord*=d; ycoord*=d; zcoord*=d; 
	}
	void normalize()
	{
		double s=this->norm();
		if (s>0.)
		{
			xcoord/=s; ycoord/=s; zcoord/=s;
		}
	}
	double norm()
	{
		return sqrt(xcoord*xcoord+ycoord*ycoord+zcoord*zcoord);
	}
	double sp(const Vector &v) const
	{
		return (xcoord*v.getx()+ycoord*v.gety()+zcoord*v.getz());
	}
	Vector vp(const Vector &v) const
	{ 
		return Vector(
			ycoord*v.getz()-zcoord*v.gety(),
			zcoord*v.getx()-xcoord*v.getz(),
			xcoord*v.gety()-ycoord*v.getx()); 
	}
};

// The class contains the value of the non-stationary signal in the receivers and I/O functions
class Picket
{
public:
	/*!   */
	struct EdsVal
	{
		double n, a, s; //!< , ,  
		EdsVal()
		{
		}
		EdsVal(const double& _n, const double& _a, const double& _s)
		{
			n=_n; 
			a=_a;
			s=_s;
		}
		EdsVal(const double& _v)
		{
			n=_v; 
			a=_v;
			s=_v;
		}
		EdsVal(const EdsVal& _e)
		{
			n=_e.n; 
			a=_e.a;
			s=_e.s;
		}
		void operator+=(const EdsVal& _v)
		{
			n+=_v.n;
			a+=_v.a;
			s+=_v.s;
		}
		EdsVal operator+(const EdsVal& _v)
		{
			return EdsVal(n+_v.n, a+_v.a, s+_v.s);
		}
		EdsVal operator-(const EdsVal& _v)
		{
			return EdsVal(n-_v.n, a-_v.a, s-_v.s);
		}
		void operator*=(const double& _m)
		{
			n*=_m;
			a*=_m;
			s*=_m;
		}
		EdsVal operator*(const double& _m)
		{
			return EdsVal(n*_m, a*_m, s*_m);
		}
		EdsVal operator/(const double& _m)
		{
			return EdsVal(n/_m, a/_m, s/_m);
		}
	};
private:
	Point3D picpnt;				//!<  
	vector<double> pictimes;	//!< 
	vector<EdsVal> piceds;		//!<  
public:
	Picket() 
	{
	}
	~Picket() 
	{
	}
	char leg[1024];
	vector<double>& GetTimes() { return pictimes; }
	vector<EdsVal>& GetEds() { return piceds; }
	Point3D& GetPnt() { return picpnt; }
	const vector<double>& GetTimes() const { return pictimes; }
	const vector<EdsVal>& GetEds() const { return piceds; }
	const Point3D& GetPnt() const { return picpnt; }
	void Clear()
	{
		pictimes.clear();
		piceds.clear();
		picpnt.x()=picpnt.y()=picpnt.z()=0.0;
	}
	/*!    edsall */
	void ReadPicket(ifstream& inf, char *buf, bool add=true, double coeff=1)
	{
		int l=0, x_int, y_int;
		double t;
		EdsVal edsv;
		{
			istringstream is(buf);
			is>>x_int>>y_int;
			if (is.fail())
				throw logic_error("bad all-file");
		}
		inf.getline(buf, 1000, '\n');
		{
			istringstream is(buf);
			is>>picpnt.x()>>picpnt.y()>>picpnt.z();
			if (is.fail())
				throw logic_error("bad all-file");
		}
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf.getline(leg, 1000, '\n');
		while (!inf.eof())
		{
			inf.getline(buf, 1000, '\n');
			if (inf.eof()) break;
			{
				istringstream is(buf);
				is>>t>>edsv.n>>edsv.a>>edsv.s;
				if (is.fail())
					break;
			}
			edsv*=coeff;
			if (add) 
			{
				pictimes.push_back(t);
				piceds.push_back(edsv);
			}
			else
				piceds[l]+=edsv;
			l++;
		}
	}
	void WritePicket(ofstream& outf)
	{
		outf<<(int)picpnt.getx()<<' '<<(int)picpnt.gety()<<endl;
		outf<<picpnt.getx()<<' '<<picpnt.gety()<<' '<<picpnt.getz()<<endl;
		outf<<"1"<<endl<<"1 1 1 1"<<endl<<"1"<<endl<<leg<<endl;
		for (int i=0; i<(int)pictimes.size(); i++)
			outf<<pictimes[i]<<" \t"<<piceds[i].n<<" \t"<<piceds[i].a<<" \t"<<piceds[i].s<<endl;
	}
	/*!      t */
	bool GetCurveValue(const double& t, EdsVal& val) const
	{
		double ht, hta, htb;
		for (int i=0; i<int(piceds.size()-1); i++)
		{
			const EdsVal& va=piceds[i];
			const EdsVal& vb=piceds[i+1];
			const double& vat=pictimes[i];
			const double& vbt=pictimes[i+1];
			if (vat<=t && t<=vbt)
			{
				ht=vbt-vat;
				hta=t-vat;
				htb=vbt-t;
				val.n=(va.n*htb+vb.n*hta)/ht;
				val.a=(va.a*htb+vb.a*hta)/ht;
				val.s=(va.s*htb+vb.s*hta)/ht;
				return true;
			}
		}
		return false;
	}
	int HevToImp2(double Pimp, double Ppause, int NT2, int ImpType)
	{
		int i,j,k,ktime,Np,Nb,ipi,ipp;
		EdsVal tmp;
		vector<int> Mask;
		vector<double> time_o;
		vector<EdsVal> val_o;

		val_o.resize(pictimes.size());

		ktime=(int)pictimes.size();
		if (ktime!=(int)pictimes.size())
			return 1;

		if(NT2<1)return 0;

		Np=2*NT2-1;
		Mask.resize(Np);
		for(k=0;k<Np;k++){
			if(ImpType){
				j=k%2;
			}
			else{
				j=(k%4)/2;
			}
			Mask[k]=2*j-1;
		}

		time_o.resize(ktime);

		for (i=0; i<ktime; i++)
			val_o[i]=piceds[i];

		Nb=ktime;
		tmp=0;

		ipi=ipp=0;
		for(k=0;k<Np;k++)
		{
			if(ipi==ipp)
				ipi++;
			else 
				ipp++;
			for(i=0;i<ktime;i++)
				time_o[i]=pictimes[i]-(ipi*Pimp+ipp*Ppause);
			j=0;
			for(i=0;i<ktime;i++)
			{
				if(pictimes[i]<=time_o[ktime-1])
				{
					while(pictimes[i]>time_o[j])
						j++;
					tmp=((piceds[j]*(pictimes[i]-time_o[j-1])+piceds[j-1]*(time_o[j]-pictimes[i]))/
						(time_o[j]-time_o[j-1]))*Mask[k];
					val_o[i]+=tmp;
				}
				else
				{
					Nb=i;
					break;
				}
			}
		}

		for(i=0;i<(int)pictimes.size();i++)
			piceds[i]=0;

		for(i=0;i<Nb;i++)
		{
			if(pictimes[i]>Ppause)
			{
				if(i>0 && Ppause/pictimes[i-1]>1.001)
				{
					double cff=(Ppause-pictimes[i-1])/(pictimes[i]-pictimes[i-1]);
					pictimes[i]=Ppause;
					piceds[i]=val_o[i]*cff+val_o[i-1]*(1.0-cff);
					i++;
				}
				break;
			}
			piceds[i]=val_o[i];
		}

		piceds.resize(i);
		pictimes.resize(i);

		return 0;
	}
};
