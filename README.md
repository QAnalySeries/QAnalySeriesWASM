# QAnalySeriesWASM
WASM / Web version of QAnalySeries tool. QAnalySeries is a partial replacement for the original AnalySeries time series analysis and correlation tool of Didier Paillard.

Created by Heiko Pälike & Sergey Kotov, last modified on 2024-03-28

We introduce a computer program for time series tuning and analysis.
The well known in paleo-climatic community program AnalySeries after Paillard et al. is restricted to the Mac OS (32-bit) and no longer compatible with current versions of Macintosh OS. QAnalySeries is an attempt to re-implement the major functionality of AnalySeries thus providing the community with a useful tool. QAnalySeries is written using Qt SDK as a free software and can be run on Macintosh, Windows and Linux systems. 
If you find this useful, please cite 
QAnalySeries-a cross-platform time series tuning and analysis tool
S Kotov, H Pälike - AGU Fall Meeting Abstracts, PP53D-1230, 2018

Data, target, time model must be a text two columns file with delimiters (space, tab, ",", ";" ). First column must be depth or time in increasing order. Adding and deleting of tie points are the same as in original AnalySeries with "Shift" and "Alt" buttons ("Shift-Alt-LeftButtonClick" to delete a tie on Ubuntu).
An experimental new version of QAnalySeries is now also available for testing at https://paloz.marum.de/QAnalySeriesWASM/index.html.
This version runs inside the user's browser, and is generated as a Webassembly "Web-App".

Limitation on mobile devices
Note that this version does not yet implement touch-screen gestures, and thus it is not currently possible to pinch-zoom. On Desktop machines a mouse-wheel allows horizontal zooming (or full zooming when additionally the Shift key is pressed).
This version includes an initial version of a Spectral Analysis function (under the Tools→Spectral Analysis Menu).
Data can be uploaded and downloaded to the local computer. Astronomical data can be downloaded on demand.
If you host this application locally, please adjust the network location to load the supplementary laskar.db3 file (included under resources), to reduce initial loading time.

This software acknowledges use of:
- the Eigen librarx (https://eigen.tuxfamily.org/index.php) MPL2.0 license
- QCustomPlot (https://www.qcustomplot.com), GPL license
- Qt WASM Toolchain (www.qt.com).

and is therefore also licensed under the GPL-3 license (see LICENSE).
