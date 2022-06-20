"""
@author : Milinda Fernando
@brief  : Compute a directed graph from a sympy expression. 
    1. Rename custom functions to symbols
    2. Graphs are merged from multiple expressions
NetworkX is required, and used to store the graph. 
"""

import sympy as sympy
import networkx as nx

class ExpressionGraph():

    def __init__(self):
        self._nx_graphs=dict()
        self._sympy_expr=dict()


    """
    Preorder traversal of the expression converting undefined functions to sympy symbols
    """
    def __pre_traversal_1(self,expr,node_list,edge_list):
        if isinstance(expr.func, sympy.core.function.UndefinedFunction):
            sym_name=str(expr.func)
            for a in expr.args:
                sym_name = sym_name + '_' + str(a)
            
            node_list.append(sympy.Symbol(sym_name))
        else:
            node_list.append(expr)
        
        for arg in expr.args:
            if isinstance(arg.func, sympy.core.function.UndefinedFunction):
                f=arg.func
                sym_name=str(f)
                for a in arg.args:
                    sym_name = sym_name + '_' + str(a)
                
                node_list.append(sympy.Symbol(sym_name))
                edge_list.append((expr,sympy.Symbol(sym_name)))
            else:
                edge_list.append((expr,arg))
                self.__pre_traversal_1(arg,node_list,edge_list)

    """
    Keep undefined function references as it is pruning but not renaming. 
    """

    def __pre_traversal_2(self,expr,node_list,edge_list):
        
        node_list.append((hash(expr), {"func": expr.func , "args": expr.args, "eval":False}))
        for arg in expr.args:
            node_list.append((hash(arg), {"func": arg.func, "args": arg.args, "eval":False}))
            edge_list.append((hash(expr),hash(arg)))
            self.__pre_traversal_2(arg,node_list,edge_list)

    """
    Pre order traversal and returns the node list and edge list
    """

    def __preorder_traversal__(self,expr):
        expr_list=list()
        edge_list=list()
        self.__pre_traversal_2(expr,expr_list,edge_list)
        return [expr_list,edge_list]

    """
    
    Generate a networkx graph for a given expression

    """
    def add_expression(self,expr,expr_name):

        G=nx.DiGraph(vname=str(expr_name))
        G.name=str(expr_name)

        # total nodes on the graph. 
        [node_list,edge_list]=self.__preorder_traversal__(expr)

        G.add_nodes_from(node_list)
        G.add_edges_from(edge_list)

        self._nx_graphs[str(expr_name)]=G
        self._sympy_expr[str(expr_name)]=expr
        
        # print("Nodes in the Graph")
        # for sub_expr in node_list:
        #     print(sub_expr)

        # print("Edges in the Graph")
        # for e in edge_list:
        #     print(e[0])
        #     print("|\n|\n")
        #     print(e[1])
        #     print("\n\n")


        #print("nodes")
        #print(G.nodes(data=True))

        #print("edges")
        #print(G.nodes(data=True))

        #nx.draw_networkx(G,pos=nx.planar_layout(G),font_size=8)
        #plt.show()

    """
    
    Adds list of sympy expressions

    """

    def add_expressions(self,outs,vnames,suffix_idx="[pp]"):
        mi = [0, 1, 2, 4, 5, 8]
        midx = ['00', '01', '02', '11', '12', '22']
        
        num_e = 0
        for i, e in enumerate(outs):
            if type(e) == list:
                num_e = num_e + len(e)
                for j, ev in enumerate(e):
                    expr_name = vnames[i] + "" + str(j) + str(suffix_idx)
                    #print("processing expr : %d var name %s[%s]" %(i,vnames[i],str(j)))
                    self.add_expression(ev,expr_name)
            elif type(e) == sympy.Matrix:
                num_e = num_e + len(e)
                for j, k in enumerate(mi):
                    expr_name = vnames[i] + "" +str(midx[j]) + str(suffix_idx)
                    #print("processing expr : %d var name %s[%s]" %(i,vnames[i],midx[j]))
                    self.add_expression(e[k],expr_name)
            else:
                num_e = num_e + 1
                #print("processing expr : %d var name %s" %(i,vnames[i]))
                expr_name = vnames[i] + str(suffix_idx)
                self.add_expression(e,expr_name)
    
    """
    compute the composed graph for all the added expressions
    """
    def composed_graph(self,verbose=False):
        g_list = list()
        for (v,g) in self._nx_graphs.items():

            if verbose:
                print("graph for var: %s" %str(v))
                print(nx.info(g))

            g_list.append(g)
        
        self._G_=nx.compose_all(g_list)
        self._G_.name="full graph added expressions"
        if verbose:
            print("full graph")
            print(nx.info(self._G_))
        
        return self._G_

    """
    plots the adjacency matrix for the composed big G
    """
    def plot_adjmatrix(self):
        import matplotlib.pyplot as plt
        A=nx.adjacency_matrix(self._G_)
        plt.spy(A,markersize=3)
        plt.show()

    """
    plots the graph for a given sympy expression
    """    
    
    def draw_graph(self,expr_name):
        import matplotlib.pyplot as plt
        g = self._nx_graphs[expr_name]
        nx.draw_networkx(g,pos=nx.planar_layout(g),font_size=6)
        plt.show()

    """
    Get the computational graph for a given expression
    """
    def get_graph(self,expr_name):
        g = self._nx_graphs[expr_name]
        return g

                








    
    