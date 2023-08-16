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
 *  This file contains the headers for formatting I/O
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
class In_Out
{
private:
	int n_of_rec;
	int len_of_rec;
	int size_of_array;

	int type_of_values;

	char input_fname[100];
	char output_fname[100];

public:
	In_Out();
	~In_Out();

	int Write_Txt_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);
	int Write_Bin_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);
	int Read_Txt_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);
	int Read_Bin_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);

	int Write_Txt_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record);
	int Write_Bin_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record);
	int Read_Txt_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record); 
	int Read_Bin_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record); 
	int Read_Bin_File_Of_Short(char *fname, short *massiv, int n_of_records, int len_of_record);

	int Read_Long_From_Txt_File(char *fname, int *number);
	int Read_Double_From_Txt_File(char *fname, double *number);
	int Write_Double_To_Txt_File(char *fname, double number);
	int Write_Long_To_Txt_File(char *fname, int number);

	int Convert_File_Of_Double_From_Bin_To_Txt(char *file_in, char *file_out);
	int Convert_File_Of_Double_From_Bin_To_Txt(char *file_in, char *file_out, int size_of_array);
	int Convert_File_Of_Double_From_Bin_To_Txt(char *file_in, char *file_out, int n_of_records, int len_of_record);

	int Convert_File_Of_Long_From_Bin_To_Txt(char *file_in, char *file_out);	
	int Convert_File_Of_Long_From_Bin_To_Txt(char *file_in, char *file_out, int size_of_array);
	int Convert_File_Of_Long_From_Bin_To_Txt(char *file_in, char *file_out, int n_of_records, int len_of_record);



	int Menu(char *input_fname, char *output_fname);



	int Write_jg(char *fname, int *ig, int *jg, int n);

    int Write_kuslau(char *fname, int n, double eps, int maxiter);
	int Write_kuslau_block(char *fname, int n_block, double eps, int maxiter);

	int Read_inftry(char *fname, int *kuzlov, int *kpar, int *l1);
	int Write_inftry(char *fname, int kuzlov, int kpar, int l1);

	int Read_1d_data(int n_1d, double *coords_1d, double *sin_1d, double *cos_1d);

	int Read_mesh_from_geoprep();

	int Read_n_eps_maxiter(char *fname, int *n, double *eps, int *maxiter);
	int Read_n_maxiter_eps(char *fname, int *n, int *maxiter, double *eps);
	int Read_inf2tr(char *fname, int *kuzlov, int *ktr, int *l1);
	int Write_inf2tr(char *fname, int kuzlov, int ktr, int l1);
};
