/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

/* IMPLEMENTE AS FUNCOES INDICADAS
 * DIGITE SEU RA: 145980
 * SUBMETA SOMENTE ESTE ARQUIVO */

#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include <gurobi_c++.h>
#include "mygraphlib.h"
#include "color_exact_heuristic.h"

int colorNaive(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit);

//------------------------------------------------------------------------------
int colorExact(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit)
/* SUBSTITUA O CONTEÚDO DESTA FUNÇÃO POR SUA IMPLEMENTAÇÃO DO ALGORITMO EXATO.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DA FUNÇÃO. */
{
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	
	GRBVar z;
	vector<GRBVar> y(gd.n);
	ListGraph::NodeMap< vector<GRBVar> > x(gd.g, vector<GRBVar>(gd.n));
	
	for (NodeIt e(gd.g); e!=INVALID; ++e) {
		for (int i = 0; i < gd.n; i++) {
			x[e][i] = model.addVar(0, 1, 0, GRB_BINARY);
		}
	}
	
	for (int i = 0; i < gd.n; i++) {
		y[i] = model.addVar(0, 1, 0, GRB_BINARY);
	}
	
	z = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
	
	model.update();
	
	for(NodeIt e(gd.g); e!=INVALID; ++e){
		GRBLinExpr expr;
		
		for (int i = 0; i < gd.n; i++) {
			expr += x[e][i];
			
			GRBLinExpr expr2;
			expr2 += x[e][i];
			model.addConstr(expr2 <= y[i]);
		}
		
		model.addConstr(expr  == 1.0);
	}
	
	for(EdgeIt e(gd.g); e!=INVALID; ++e){
		for (int i = 0; i < gd.n; i++) {
			GRBLinExpr expr;
			
			expr += x[gd.g.u(e)][i];
			expr += x[gd.g.v(e)][i];
			
			model.addConstr(expr  <= 1.0);
		}
	}
	
	GRBLinExpr expr;
	for (int i = 0; i < gd.n; i++) {
		expr += y[i];
	}
	expr = expr - z;
	model.addConstr(expr == 0);
	
	model.update();

	GRBLinExpr exprObj = z;
	model.setObjective(exprObj, GRB_MINIMIZE);
	model.update();
	
	if (timeLimit > 0){
		model.getEnv().set(GRB_DoubleParam_TimeLimit,timeLimit);
	}
	
	try {
		model.optimize();
		int status = model.get(GRB_IntAttr_Status);
		
		for (int i = 0; i < gd.n; i++) {
			for(NodeIt e(gd.g); e!=INVALID; ++e){
				if( x[e][i].get(GRB_DoubleAttr_X) == 1.0){
					color[e] = i+1;
				}
			}
		}
		
		if (status == GRB_OPTIMAL){
			lowerBound = z.get(GRB_DoubleAttr_X);
			upperBound = z.get(GRB_DoubleAttr_X);
			return 1;
		}else{
			lowerBound = ceil(model.get(GRB_DoubleAttr_MinBound));
			upperBound = z.get(GRB_DoubleAttr_X);
			return 0;
		}
	}catch(GRBException e) {
		cerr << "Nao foi possivel resolver o PLI." << endl;
		cerr << "Codigo de erro = " << e.getErrorCode() << endl;
		cerr << e.getMessage();
	}
	
	return 0;
}
//------------------------------------------------------------------------------
int colorHeuristic(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit)
/* SUBSTITUA O CONTEÚDO DESTA FUNÇÃO POR SUA IMPLEMENTAÇÃO DA HEURISTICA.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DA FUNÇÃO. */
{
	clock_t begin = clock();
	
	bool inteiro = false;
	Node lastNode;
	int lastColor = -1;
	int lastConst = -1;
	GRBConstr lastConstr;
	
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	
	GRBVar z;
	vector<GRBVar> y(gd.n);
	ListGraph::NodeMap< vector<GRBVar> > x(gd.g, vector<GRBVar>(gd.n));
	
	for (NodeIt e(gd.g); e!=INVALID; ++e) {
		for (int i = 0; i < gd.n; i++) {
			x[e][i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
	}
	
	for (int i = 0; i < gd.n; i++) {
		y[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}
	
	z = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
	
	model.update();
	
	for(NodeIt e(gd.g); e!=INVALID; ++e){
		GRBLinExpr expr;
		for (int i = 0; i < gd.n; i++) {
			expr += x[e][i];
			
			GRBLinExpr expr2;
			expr2 += x[e][i];
			model.addConstr(expr2 <= y[i]);
		}
		model.addConstr(expr  == 1.0);
	}
	
	for(EdgeIt e(gd.g); e!=INVALID; ++e){
		for (int i = 0; i < gd.n; i++) {
			GRBLinExpr expr;
			
			expr += x[gd.g.u(e)][i];
			expr += x[gd.g.v(e)][i];
			
			model.addConstr(expr  <= 1.0);
		}
	}
	
	GRBLinExpr expr;
	for (int i = 0; i < gd.n; i++) {
		expr += y[i];
	}
	expr = expr - z;
	model.addConstr(expr == 0);
	
	model.update();
	
	if (timeLimit > 0){
		model.getEnv().set(GRB_DoubleParam_TimeLimit,timeLimit);
	}
	
	do{
		try {
			GRBLinExpr exprObj = z;
			model.setObjective(exprObj, GRB_MINIMIZE);
			model.optimize();
			int status = model.get(GRB_IntAttr_Status);
			
			clock_t posOpt = clock();
	
			if( (double(posOpt - begin) / CLOCKS_PER_SEC) > timeLimit ){
				for (int i = 0; i < gd.n; i++) {
					for(NodeIt e(gd.g); e!=INVALID; ++e){
						if( x[e][i].get(GRB_DoubleAttr_X) == 1.0){
							color[e] = i+1;
						}
					}
				}
				lowerBound = model.get(GRB_DoubleAttr_MinBound);
				upperBound = z.get(GRB_DoubleAttr_X);
				return 0;
			}
			
			if (timeLimit > 0){
				model.getEnv().set(GRB_DoubleParam_TimeLimit, timeLimit - ((posOpt - begin) / CLOCKS_PER_SEC));
			}
			
			if(status != GRB_OPTIMAL){
				GRBLinExpr constrAdd = x[lastNode][lastColor];
				if(lastConst == 1){
					model.remove(lastConstr);
					model.update();
					lastConstr = model.addConstr(constrAdd == 0.0);
					lastConst = 0;
				}else{
					model.remove(lastConstr);
					model.update();
					lastConstr = model.addConstr(constrAdd == 1.0);
					lastConst = 1;
				}
				model.update();
			}else{
				inteiro = true;
				for (int i = 0; i < gd.n; i++) {
					for(NodeIt e(gd.g); e!=INVALID; ++e){
						double valor = x[e][i].get(GRB_DoubleAttr_X);
						if( valor != 1.0 && valor != 0.0){
							inteiro = false;
							if(valor < 0.5){
								GRBLinExpr constrAdd = x[e][i];
								lastConstr = model.addConstr(constrAdd == 0.0);
								lastNode = e;
								lastColor = i;
								lastConst = 0;
								model.update();
							}else{
								GRBLinExpr constrAdd = x[e][i];
								lastConstr = model.addConstr(constrAdd == 1.0);
								lastNode = e;
								lastColor = i;
								lastConst = 1;
								model.update();
							}
						}
						break;
					}
				}
				
				if(inteiro){
					for (int i = 0; i < gd.n; i++) {
						for(NodeIt e(gd.g); e!=INVALID; ++e){
							if( x[e][i].get(GRB_DoubleAttr_X) == 1.0){
								color[e] = i+1;
							}
						}
					}
					lowerBound = z.get(GRB_DoubleAttr_X);
					upperBound = z.get(GRB_DoubleAttr_X);
					return 1;
				}
			}
		}catch(GRBException e) {
			GRBLinExpr constrAdd = x[lastNode][lastColor];
			if(lastConst == 1){
				model.remove(lastConstr);
				model.update();
				lastConstr = model.addConstr(constrAdd == 0.0);
				lastConst = 0;
			}else{
				model.remove(lastConstr);
				model.update();
				lastConstr = model.addConstr(constrAdd == 1.0);
				lastConst = 1;
			}
			model.update();
		}
	}while(!inteiro);
	
	return colorNaive(gd, color, lowerBound, upperBound, timeLimit);
}
//------------------------------------------------------------------------------
int colorNaive(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit)
/* Algoritmo ingênuo para o problema de coloração de vértices */
{
   lowerBound = 1;
   int next = lowerBound;
	
	for(NodeIt i(gd.g); i != INVALID; ++i){
      color[i] = next;
      next++;
	}
	next--;
	upperBound = next;
	
	return 0;
}
//------------------------------------------------------------------------------

