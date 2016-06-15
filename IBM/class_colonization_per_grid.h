#ifndef __COLONIZATION_H_
#define __COLONIZATION_H_

class ColonizationPerGrid{
	private:
		int ngrid;
		double *grid_col;

	public:
		ColonizationPerGrid(int _ngrid){
			this->ngrid=_ngrid;
			this->grid_col = new double[_ngrid * _ngrid];
			for ( int i = 0 ; i< (this->ngrid)*(this->ngrid) ; i++){
				this->grid_col[i]=0.0;
			}
		}

		~ColonizationPerGrid(){
			delete[] this->grid_col;
		}

		int n_cell(){
			return this->ngrid;
		}

		double get_grid(int i, int j){
			return this->grid_col[i + j*ngrid];
		}

		/*double get_grid_direct(int i){
			return this->grid_col[i];
		}*/

		vector<int> rand_choose_gridC(double totalcol){
			//assert(totalcol>0);
			//double rand = get_random();
			double temptotal = totalcol*get_random();
			//assert(temptotal<totalcol);
			double sum = 0.0;
			vector<int> result(2);

			/*for (int i = 0 ; i < ngrid ; i++){
				for (int j = 0 ; j < ngrid ; j++){
					//cout<<this->grid_col[i + j*ngrid]<<endl;
					assert(0.0<= this->grid_col[i + j*ngrid]);
					sum+= this->grid_col[i + j*ngrid];
					//if( temptotal <= sum ){
					//	result[0] = i;
					//	result[1] = j;
					//	goto end;
					//}
				}
			}*/
			//if(sum!= totalcol){cout<<"sum: "<<sum<<"totalcol: "<<totalcol - sum;};
			//assert(sum == totalcol);
			//sum = 0.0;
			for (int i = 0 ; i < ngrid ; i++){
				for (int j = 0 ; j < ngrid ; j++){
					//cout<<this->grid_col[i + j*ngrid]<<endl;
					assert(0.0<= this->grid_col[i + j*ngrid]);
					sum+= this->grid_col[i + j*ngrid];
					if( temptotal <= sum ){
						result[0] = i;
						result[1] = j;
						goto end;
					}
				}
			}
			end:
			//if(grid_col[result[0] + result[1]*ngrid] == 0){cout<<"rand:"<< rand <<"tempcol: " << temptotal<<", totalcol: "<<totalcol;};
			//assert(grid_col[result[0] + result[1]*ngrid] != 0);
			return result;
		}

		void modify_col_rate(int i, int j, double col){
			this->grid_col[i+j*ngrid]+= col;
			if(this->grid_col[i+j*ngrid]<0.0 || this->grid_col[i+j*ngrid]<chop_0){
				this->grid_col[i+j*ngrid] = 0.0;
			}
		}
};

ColonizationPerGrid *initializeColonizationGrid(Grid *g){
	int ncell= g->get_number_of_cells();
	double col_rate=0.0;

	ColonizationPerGrid *Cgrid = new ColonizationPerGrid(ncell);

	for( int i = 0 ; i < ncell ; i++){
		for( int j = 0 ; j<ncell ; j++ ){
			col_rate = g->colrate_per_grid(i,j);
			Cgrid->modify_col_rate(i,j,col_rate);
		}
	}
	return Cgrid;
}

void outputcpg(ColonizationPerGrid *cpg){
	double temp_col, tot_col=0.0;
	for( int i=0; i<cpg->n_cell(); i++){
		for( int j=0; j<cpg->n_cell(); j++){
			temp_col = cpg->get_grid(i,j);
			tot_col+= temp_col;
			cout<<" in grid " << i<< " and "<<j<<" the col rate is "<< temp_col<<"\n";
		}
	}
	cout<<" total colonization cpg  "<< tot_col<<endl;
}

			/*int chosen_grid, ni, nj;
			for( int i =0; i<20; i++){
				chosen_grid = cpg->rand_choose_gridC(total_colonization);
				ni= (int) chosen_grid%ncell;
				nj= (int)(chosen_grid - ni)/ncell;
				cout<<" the chosen cell is " << ni<<" , "<< nj<<endl;
			}*/

#endif // __COLONIZATION_H_