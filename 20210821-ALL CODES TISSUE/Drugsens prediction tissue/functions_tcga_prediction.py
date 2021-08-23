# Last updated on Aug 21th, 2021,  fjq19@mails.tsinghua.edu.cn
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
import plotly.graph_objects as go
import plotly.io as pio
import random

def my_pca_RNA_2(X1,X2,y,y_map, pca_dim=2):
    # pca
    pca=PCA(n_components=min(pca_dim, X1.shape[0]))
    X_pca_1 = pca.fit_transform(X1)
    print('sum(pca.explained_variance_ratio_)',sum(pca.explained_variance_ratio_))
    print(X_pca_1.shape)
    X_pca_2 = pca.transform(X2)
    print('sum(pca.explained_variance_ratio_)',sum(pca.explained_variance_ratio_))
    print(X_pca_2.shape)
    
    X_pca = np.r_[X_pca_1,X_pca_2]

    plt.scatter(X_pca_1[:,0], X_pca_1[:,1],c=y_map[0:len(X1)], cmap='jet',s=y[0:len(X1)]+20)
    plt.colorbar()
    plt.scatter(X_pca_2[:,0], X_pca_2[:,1],c=y_map[len(X1):len(X1)+len(X2)], cmap='gnuplot',s=y[len(X1):len(X1)+len(X2)]+20,marker = '*')
    plt.colorbar()
    plt.show()
    return X_pca_1,X_pca_2

def Classification_By_Knn(X_train, y_train, X_test, neigh_num=3):
    neigh = KNeighborsClassifier(n_neighbors=neigh_num,weights='distance')
    neigh.fit(X_train, y_train)
    y_test=neigh.predict(X_test)
    y_test_proba = neigh.predict_proba(X_test)
    return y_test,y_test_proba

## Probability Estimates with go.Contour
def KNN_prob_estimate_plot_3(X,y,X_test,y_test,y_test_ic50,neighbor_k,thresh=100,save_path=''):
    mesh_size = .02
    margin = 0.5
    
    # Create a mesh grid on which we will run our model
    x_min, x_max = X[:, 0].min() - margin, X[:, 0].max() + margin
    y_min, y_max = X[:, 1].min() - margin, X[:, 1].max() + margin
    xrange = np.arange(x_min, x_max, mesh_size)
    yrange = np.arange(y_min, y_max, mesh_size)
    xx, yy = np.meshgrid(xrange, yrange)

    # Create classifier, run predictions on grid
    clf = KNeighborsClassifier(neighbor_k, weights='distance')
    clf.fit(X, y)
    Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]
    Z = Z.reshape(xx.shape)
    Z_max=np.max(Z)
    Z[np.where(Z>Z_max*0.7)]=Z_max*0.7
    Z_min=np.min(Z)
    Z[np.where(Z<Z_min*0.7)]=Z_min*0.7

    layout = go.Layout(
#         autosize=True,
        autosize=False,
        width=1000,
        height=800
    )
    trace_specs = [
        [X, y, 0, 'Train', 'square'],
        [X, y, 1, 'Train', 'circle']
    ]

    fig = go.Figure(data=[
        go.Scatter(
            x=X[y==label, 0], y=X[y==label, 1],
            name=f'{split} Split, Label {label}',
            mode='markers', marker_symbol=marker
        )
        for X, y, label, split, marker in trace_specs
    ],layout = layout)
    fig.update_traces(
        marker_size=14, marker_line_width=1.5,
        marker_color="lightyellow"
    )
    
    fig.add_trace(
        go.Scatter(
            x=X_test[y_test_ic50>thresh, 0], y=X_test[y_test_ic50>thresh, 1],name='Test-IC50>100',
            mode='markers', marker_symbol='star',
            marker_size=25, marker_line_width=1.5,marker_color="red",marker_line_color='red'
        )
    )
    fig.add_trace(
        go.Scatter(
            x=X_test[y_test_ic50<thresh, 0], y=X_test[y_test_ic50<thresh, 1],name='Test-IC50<100',
            mode='markers', marker_symbol='star',
            marker_size=25, marker_line_width=1.5,marker_color="darkblue",marker_line_color='darkblue'
        )
    )
    fig.add_trace(
        go.Heatmap(
            x=xrange,
            y=yrange,
            z=-Z,
            showscale=True,
            colorscale='rdylgn',
            opacity=1,
            name='Score',
            type='heatmap'
        )
    )
    fig.update_xaxes(color='black',visible=True,title='PC 1',ticks='outside',tickcolor='black',tickfont_size=22,tickfont_family='Arial',tickwidth=2,showline=True,linecolor='black',linewidth=2)
    fig.update_yaxes(color='black',visible=True,title='PC 2',ticks='outside',tickcolor='black',tickfont_size=22,tickfont_family='Arial',tickwidth=2,dtick=1,showline=True,linecolor='black',linewidth=2)
    if len(save_path)>0:
        pio.write_image(fig,save_path+ '.png')
        pio.write_image(fig,save_path+ '.svg')
    fig.show()
    return 0

