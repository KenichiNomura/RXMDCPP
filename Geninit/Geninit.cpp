#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <assert.h>
#include "Geninit.h"

Geninit::Geninit(char *fname, int repX, int repY, int repZ, int procX, int procY, int procZ)
{

	pi = atan(1.0)*4.0;

	mc[0] = repX; mc[1] = repY; mc[2] = repZ; 
	mctot = mc[0]*mc[1]*mc[2];
	vprocs[0] = procX; vprocs[1] = procY; vprocs[2] = procZ; 
	nprocs = vprocs[0]*vprocs[1]*vprocs[2];

	std::ifstream inputFileContainer;
	inputFileContainer.open(fname);

	int lineNum = 0;
	std::string tempVarToReadLine;
	if (inputFileContainer.is_open())
	{
		while (getline(inputFileContainer, tempVarToReadLine))
		{
			std::istringstream iss(tempVarToReadLine);
			if (lineNum == 0)
			{
				iss >> natoms >> fnote;
				ctype0 = new char*[natoms];
				ctype1 = new char*[natoms*mctot];

				lnatoms  = new int[nprocs];
				lnatoms1 = new int[nprocs];
				lnatoms2 = new int[nprocs];
				
				itype0 = new int[natoms];
				itype1 = new double[natoms*mctot];
			
				pos0 = new double*[natoms];
				pos1 = new double*[natoms*mctot];
				rr3 =  new double*[natoms*mctot];
				
				for (int i = 0; i < natoms; i++)  
				{     
					pos0[i] = new double[3];
					ctype0[i] = new char[256];
				}
				for (int i = 0; i < natoms*mctot; i++) 
				{
					pos1[i] = new double[3];
					rr3[i]  = new double[3];
					ctype1[i] = new char[256];
				}				
			}
			else if (lineNum == 1)
				iss >> L1 >> L2 >> L3 >> Lalpha >> Lbeta >> Lgamma;
			else
			{
				iss >> ctype0[lineNum-2] >> pos0[0][lineNum-2] >> pos0[1][lineNum-2] >> pos0[2][lineNum-2];
				if (ctype0[lineNum-2] == "C") itype0[lineNum-2] = 1;
				else if (ctype0[lineNum-2] == "H") itype0[lineNum-2] = 2;
				else if (ctype0[lineNum-2] == "O") itype0[lineNum-2] = 3;
				else if (ctype0[lineNum-2] == "N") itype0[lineNum-2] = 4;
			}
		
			lineNum++;
		}
	}
}

void Geninit::HMatrix()
{
	double lal, lbe, lga, hh1, hh2;

	lal = Lalpha*pi/180.0; lbe = Lbeta*pi/180.0; lga = Lgamma*pi/180.0;

	hh1 = L3*(cos(lal) - cos(lbe)*cos(lga))/sin(lga);
	hh2 = L3*sqrt(1.0 - pow(cos(lal),2) - pow(cos(lbe),2) - pow(cos(lga),2) + 2.0*cos(lal)*cos(lbe)*cos(lga))/sin(lga);

	H[0][0] = L1; 		H[0][1] = L2*cos(lga); 		H[0][2] = L3*cos(lbe);
	H[1][0] = 0.0;		H[1][1] = L2*sin(lga);		H[1][2] = hh1;
	H[2][0] = 0.0;		H[2][1] = 0.0;			H[2][2] = hh2;
	
}

void Geninit::Hinv()
{
	double detm;

	Hi[0][0] = H[1][1]*H[2][2] - H[1][2]*H[2][1];
	Hi[0][1] = H[0][2]*H[2][1] - H[0][1]*H[2][2];
	Hi[0][2] = H[0][1]*H[1][2] - H[0][2]*H[1][1];
	Hi[1][0] = H[1][2]*H[2][0] - H[1][0]*H[2][2];
	Hi[1][1] = H[0][0]*H[2][2] - H[0][2]*H[2][0];
	Hi[1][2] = H[0][2]*H[1][0] - H[0][0]*H[1][2];
	Hi[2][0] = H[1][0]*H[2][1] - H[1][1]*H[2][0];
	Hi[2][1] = H[0][1]*H[2][0] - H[0][0]*H[2][1];
	Hi[2][2] = H[0][0]*H[1][1] - H[0][1]*H[1][0];

	detm = H[0][0]*H[1][1]*H[2][2] + H[0][1]*H[1][2]*H[2][0]
	     + H[0][2]*H[1][0]*H[2][1] - H[0][2]*H[1][1]*H[2][0]
	     - H[0][1]*H[1][0]*H[2][2] - H[0][0]*H[1][2]*H[2][1];
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) Hi[i][j] = Hi[i][j]/detm;

}

void Geninit::Replicate(int NormalOrReal)
{
	int ntot;
	for (int ix = 0; ix < mc[0]; ix++)
		for (int iy = 0; iy < mc[1]; iy++)
			for (int iz = 0; iz < mc[2]; iz++)
				for (int i = 0; i < natoms; i++)
				{
					if (NormalOrReal == 1)
					{
						rr[0] = Hi[0][0]*pos0[i][0] + Hi[0][1]*pos0[i][1] + Hi[0][2]*pos0[i][2];
						rr[1] = Hi[1][0]*pos0[i][0] + Hi[1][1]*pos0[i][1] + Hi[1][2]*pos0[i][2];
						rr[2] = Hi[2][0]*pos0[i][0] + Hi[2][1]*pos0[i][1] + Hi[2][2]*pos0[i][2];
					}
					rr[0] = (rr[0] + ix)/mc[0];
					rr[1] = (rr[1] + iy)/mc[1];
					rr[2] = (rr[2] + iz)/mc[2];

					pos1[ntot][0] = rr[0]; pos1[ntot][1] = rr[1]; pos1[ntot][2] = rr[2];
					ctype1[ntot] = ctype0[i];
					itype1[ntot] = itype0[i] + ntot*pow(10,-13);
					ntot++;					
				}
			
	natoms = natoms*mctot;
	L1 *= mc[1]; L2*= mc[2]; L3 *= mc[3];
	HMatrix();
	Hinv();
}

void Geninit::IntroduceVoid(double spaceX, double spaceY, double spaceZ)
{
	for (int i = 0; i < natoms; i++)
	{
		rr3[i][0] = spaceX + H[0][0]*pos1[i][0] + H[0][1]*pos1[i][1] + H[0][2]*pos1[i][2];
		rr3[i][1] = spaceY + H[1][0]*pos1[i][0] + H[1][1]*pos1[i][1] + H[1][2]*pos1[i][2];
		rr3[i][2] = spaceZ + H[2][0]*pos1[i][0] + H[2][1]*pos1[i][1] + H[2][2]*pos1[i][2];
	}

	L1 += spaceX;
	L2 += spaceY;
	L3 += spaceZ;

	HMatrix();
	Hinv();

	for (int i = 0; i < natoms; i++)
	{
		pos1[i][0] = Hi[0][0]*rr3[i][0] + Hi[0][1]*rr3[i][1] + Hi[0][2]*rr3[i][2];
		pos1[i][1] = Hi[1][0]*rr3[i][0] + Hi[1][1]*rr3[i][1] + Hi[1][2]*rr3[i][2];
		pos1[i][2] = Hi[2][0]*rr3[i][0] + Hi[2][1]*rr3[i][1] + Hi[2][2]*rr3[i][2];
	}
}

void Geninit::WrapBack(int spaceInX, int spaceInY, int spaceInZ)
{
	double minX = pow(10,7), minY = pow(10,7), minZ = pow(10,7);
	double maxX = -10000000.0, maxY = -10000000.0, maxZ = -10000000.0;

	for (int i = 0; i < natoms; i++)
	{
		if ( !spaceInX && (pos1[i][0] <= minX) ) minX = pos1[i][0];
		if ( !spaceInY && (pos1[i][1] <= minY) ) minY = pos1[i][1];
		if ( !spaceInZ && (pos1[i][2] <= minZ) ) minZ = pos1[i][2];
	}

	for (int i = 0; i < natoms; i++)
	{
		if (!spaceInX) pos1[i][0] -= minX;
		if (!spaceInY) pos1[i][1] -= minY;
		if (!spaceInZ) pos1[i][2] -= minZ;
	}

	for (int i = 0; i < natoms; i++)
	{
		if (!spaceInX) pos1[i][0] = std::fmod(pos1[i][0],1.0);
		if (!spaceInY) pos1[i][1] = std::fmod(pos1[i][1],1.0);
		if (!spaceInZ) pos1[i][2] = std::fmod(pos1[i][2],1.0);

		pos1[i][0] += pow(10,-9);
		pos1[i][1] += pow(10,-9);
		pos1[i][2] += pow(10,-9);
	}

	minX = pow(10,7); minY = pow(10,7); minZ = pow(10,7); 
	for (int i = 0; i < natoms; i++)
	{
		if (pos1[i][0] < minX) minX = pos1[i][0];
		if (pos1[i][1] < minY) minY = pos1[i][1];
		if (pos1[i][2] < minZ) minZ = pos1[i][2];

		if (pos1[i][0] > maxX) maxX = pos1[i][0];
                if (pos1[i][1] > maxY) maxY = pos1[i][1];
                if (pos1[i][2] > maxZ) maxZ = pos1[i][2];
	}

	rmin[0] = minX; rmin[1] = minY; rmin[2] = minZ;
	rmax[0] = maxX; rmax[1] = maxY; rmax[2] = maxZ;
}

void Geninit::ErrorCheck()
{

	double x,y,z, rr[3], rr1[3], dtype;
	int sID, ii;
	int idx, idy, idz;

	assert(rmin[0] > 0 && rmax[0] < 1 && "Coordinates in X direction out of bounds");
	assert(rmin[1] > 0 && rmax[1] < 1 && "Coordinates in Y direction out of bounds");
	assert(rmin[2] > 0 && rmax[2] < 1 && "Coordinates in Z direction out of bounds");

	for (int i = 0; i < nprocs; i++)
	{
		lnatoms[i]  = 0;
		lnatoms1[i] = 0;
		lnatoms2[i] = 0;
	}		

	for (int i = 0; i < natoms; i++)
	{
		x = pos1[i][0]*vprocs[0];
		y = pos1[i][1]*vprocs[1];
		z = pos1[i][2]*vprocs[2];

		sID = x + y*vprocs[0] + z*vprocs[0]*vprocs[1];
		assert(sID < nprocs && "Processor ID out of bounds");
		lnatoms[sID]++;
	}
}

void Geninit::OutputBin()
{

	 double x,y,z, rr[3], rr1[3], dtype;
        int sID, ii;
        int idx, idy, idz;

	// Count the cumulative distribution of atoms across processors
	lnatoms1[0] = 0;
	for (int i = 1; i < nprocs; i++) lnatoms1[i] = lnatoms1[i-1] + lnatoms[i-1];

	// Find reduced box Length
	lbox[0] = 1.0/static_cast<double>(vprocs[0]);
	lbox[1] = 1.0/static_cast<double>(vprocs[1]);
	lbox[2] = 1.0/static_cast<double>(vprocs[2]);

	for (int i = 0; i < nprocs; i++) lnatoms2[i] = 0.0;


	std::fstream all_binary("all.bin", std::ios::out | std::ios::in | std::ios::binary);
	// Write atoms to respective positions in all.bin binary file
	for (int i = 0; i < natoms; i++)
	{
		x = pos1[i][0]*vprocs[0];
		y = pos1[i][1]*vprocs[1];
		z = pos1[i][2]*vprocs[2];

		sID = x + y*vprocs[0] + z*vprocs[0]*vprocs[1];
		ii = lnatoms1[sID] + lnatoms2[sID];

		all_binary.seekg(ii*32, std::ios::beg);
		all_binary.write(reinterpret_cast<char*>(&pos1[i][0]), sizeof(double));
		all_binary.write(reinterpret_cast<char*>(&pos1[i][1]), sizeof(double));	
		all_binary.write(reinterpret_cast<char*>(&pos1[i][2]), sizeof(double));
		all_binary.write(reinterpret_cast<char*>(&itype1[i]), sizeof(double));
		lnatoms2[sID]++;
	}	

	all_binary.close();
	
	for (int i = 0; i < nprocs; i++)
		assert (lnatoms[i] == lnatoms2[i] && "Error on distributing atoms on different cores");

	qq = 0; vv[0] = 0; vv[1] = 0; vv[2] = 0; qfsp = 0; qfsv = 0;

	std::fstream xyz_out("geninit.xyz", std::ios::out);
	std::fstream rxff("rxff.bin", std::ios::out | std::ios::binary);
	std::fstream readBin("all.bin", std::ios::in | std::ios::binary);

	rxff.write(reinterpret_cast<char*>(nprocs), sizeof(int));
	rxff.write(reinterpret_cast<char*>(vprocs[0]),sizeof(int));
	rxff.write(reinterpret_cast<char*>(vprocs[1]),sizeof(int));
	rxff.write(reinterpret_cast<char*>(vprocs[2]),sizeof(int));
	rxff.write(reinterpret_cast<char*>(&L1), sizeof(double));
	rxff.write(reinterpret_cast<char*>(&L2), sizeof(double));
	rxff.write(reinterpret_cast<char*>(&L3), sizeof(double));
	rxff.write(reinterpret_cast<char*>(&Lalpha), sizeof(double));
	rxff.write(reinterpret_cast<char*>(&Lbeta), sizeof(double));
	rxff.write(reinterpret_cast<char*>(&Lgamma), sizeof(double));

	xyz_out << lnatoms[nprocs-1] << std::endl;
	xyz_out << std::setw(20) << L1 <<  std::setw(20) << L2 <<  std::setw(20) << L3
	        <<  std::setw(20)<< Lalpha <<  std::setw(20) << Lbeta <<  std::setw(20) << Lgamma << std::endl;

	for (int i = 0; i < nprocs; i++)
	{
		rxff.write(reinterpret_cast<char*>(&lnatoms[i]), sizeof(int));

		idx = i%vprocs[0];
		idy = (i/vprocs[0])%vprocs[1];
		idz = i/(vprocs[0]*vprocs[1]);

		obox[0] = lbox[0]*idx;
		obox[1] = lbox[1]*idy;
		obox[2] = lbox[2]*idz;

		rxff.write(reinterpret_cast<char*>(&obox[0]), sizeof(double));
		rxff.write(reinterpret_cast<char*>(&obox[1]), sizeof(double));
		rxff.write(reinterpret_cast<char*>(&obox[2]), sizeof(double));

		for (int j = 0; j < lnatoms[i]; j++)
		{
			readBin.read(reinterpret_cast<char*>(&rr[0]), sizeof(double));
			readBin.read(reinterpret_cast<char*>(&rr[1]), sizeof(double));
			readBin.read(reinterpret_cast<char*>(&rr[2]), sizeof(double));
			readBin.read(reinterpret_cast<char*>(&dtype), sizeof(double));

			rxff.write(reinterpret_cast<char*>(&rr[0]), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&rr[1]), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&rr[2]), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&vv[0]), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&vv[1]), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&vv[2]), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&qq), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&dtype), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&qfsp), sizeof(double));
			rxff.write(reinterpret_cast<char*>(&qfsv), sizeof(double));

			rr1[0] = H[0][0]*rr[0] + H[0][1]*rr[1] + H[0][2]*rr[2];
			rr1[1] = H[1][0]*rr[0] + H[1][1]*rr[1] + H[1][2]*rr[2];
			rr1[2] = H[2][0]*rr[0] + H[2][1]*rr[1] + H[2][2]*rr[2];

			int type = dtype;
			
			xyz_out << std::setw(20) << type << std::setw(20) << rr1[0] << std::setw(20) << rr1[1] << std::setw(20) << rr1[2] << std::endl;
		}
	}
	xyz_out.close();
	rxff.close();
	readBin.close();

}
