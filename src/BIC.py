import numpy as np, pandas as pd
import itertools

# from scipy import linalg
# import matplotlib.pyplot as plt
# import matplotlib as mpl

from sklearn import mixture

print(__doc__)

# load samples

X = pd.read_csv('../matrix/ten_cell_types_lncRNA.xls', sep='\t', index_col=0)

columns =[x for x in X.columns if x.find('single') == -1]

X = X[columns]

X['sum'] = X.sum(axis=1)

X = X[X['sum'] > 500]
del X['sum']
name = X['name']
del X['name']

print X.shape

lowest_bic = np.infty

bic = []
n_components_range = range(1, 11)
cv_types = ['spherical', 'tied', 'diag', 'full']

results = {}

for cv_type in cv_types:
    results[cv_type] = []
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        try:
            gmm = mixture.GaussianMixture(n_components=n_components,
                                          covariance_type=cv_type)
            gmm.fit(X)
            cur_bic = gmm.bic(X)
            bic.append(cur_bic)
            results[cv_type].append(cur_bic)
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
            print cv_type, n_components, ' Done'
        except:
            bic.append(None)
            results[cv_type].append(None)
bic = np.array(bic)
color_iter = itertools.cycle(['navy', 'turquoise', 'cornflowerblue',
                              'darkorange'])
clf = best_gmm
bars = []

# Plot the BIC scores
# if False:
#     spl = plt.subplot(2, 1, 1)
#     for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
#         xpos = np.array(n_components_range) + .2 * (i - 2)
#         bars.append(plt.bar(xpos, bic[i * len(n_components_range):
#                                       (i + 1) * len(n_components_range)],
#                             width=.2, color=color))
#     plt.xticks(n_components_range)
#     plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
#     plt.title('BIC score per model')
#     # xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
#     #     .2 * np.floor(bic.argmin() / len(n_components_range))
#     # plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
#     spl.set_xlabel('Number of components')
#     spl.legend([b[0] for b in bars], cv_types)

df = pd.DataFrame(results)
df.index = n_components_range
# df.columns = cv_types

df.to_csv('BIC.xls', sep='\t')

# Plot the winner
# splot = plt.subplot(2, 1, 2)
# Y_ = clf.predict(X)
# for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
#                                            color_iter)):
#     v, w = linalg.eigh(cov)
#     if not np.any(Y_ == i):
#         continue
#     plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)
#
#     # Plot an ellipse to show the Gaussian component
#     angle = np.arctan2(w[0][1], w[0][0])
#     angle = 180. * angle / np.pi  # convert to degrees
#     v = 2. * np.sqrt(2.) * np.sqrt(v)
#     ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
#     ell.set_clip_box(splot.bbox)
#     ell.set_alpha(.5)
#     splot.add_artist(ell)
#
# plt.xticks(())
# plt.yticks(())
# plt.title('Selected GMM: full model, 2 components')
# plt.subplots_adjust(hspace=.35, bottom=.02)
# plt.show()