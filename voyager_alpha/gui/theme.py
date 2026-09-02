APP_STYLE = r"""
* {
    font-family: "Segoe UI";
    font-size: 12px;
    color: #dbe2e5;
    letter-spacing: 0px;
}
QMainWindow, QWidget#appRoot, QWidget#workspaceRoot {
    background: #0b0f12;
}
QFrame#titleBar {
    background: #111619;
    border-bottom: 1px solid #283238;
}
QFrame#logoPlate {
    background: #080c0f; border: 1px solid #3b4a51; border-radius: 6px;
}
QLabel#brandName { font-size: 19px; font-weight: 650; color: #f2f5f6; }
QLabel#brandSub { color: #7f9099; font-size: 10px; }
QPushButton#moduleTab {
    min-width: 172px; min-height: 38px;
    background: #161d21; border: 1px solid #273138; border-bottom: 2px solid transparent;
    padding: 0 14px; color: #aebbc1;
}
QPushButton#moduleTab:hover { background: #1b2429; color: #f3f6f7; }
QPushButton#moduleTab:checked {
    background: #20292e; color: #ffffff; border-bottom: 2px solid #2bc7cf;
}
QPushButton#windowButton {
    min-width: 34px; max-width: 34px; min-height: 30px;
    background: transparent; border: 1px solid transparent; color: #a7b3b9;
}
QPushButton#windowButton:hover { background: #263036; color: #ffffff; }
QPushButton#closeButton { min-width: 38px; min-height: 30px; background: transparent; border: 0; }
QPushButton#closeButton:hover { background: #b8423c; color: white; }
QFrame#toolbar {
    background: #0e1417; border-bottom: 1px solid #263138;
}
QFrame#statusBar {
    background: #0d1316; border-top: 1px solid #334148;
}
QLabel#statusPercent {
    color: #eef9fa; background: #1d2a2f; border: 1px solid #3a4b52;
    border-radius: 3px; padding: 3px 4px; font-size: 10px; font-weight: 650;
}
QFrame#leftPanel, QFrame#rightPanel {
    background: #10171b; border: 1px solid #2a353b;
}
QFrame#centerPanel { background: #0b1013; }
QStackedWidget { background: #0b0f12; }
QFrame#section {
    background: #121a1e; border: 1px solid #2d3940; border-radius: 4px;
}
QFrame#sectionAccent {
    background: #111b1f; border: 1px solid #2f5960; border-radius: 4px;
}
QLabel#panelTitle { font-size: 14px; font-weight: 650; color: #eef3f5; }
QLabel#sectionTitle { font-size: 11px; font-weight: 650; color: #b9c7cc; }
QLabel#caption { font-size: 9px; color: #78909a; text-transform: uppercase; }
QLabel#muted { color: #80919a; }
QLabel#value { color: #edf2f4; font-weight: 550; }
QLabel#cyanValue { color: #2dd7df; font-size: 20px; font-weight: 500; }
QLabel#ok { color: #4dc56a; }
QLabel#warn { color: #efbc3f; }
QLabel#error { color: #ef6257; }
QLabel#badgeKnown { background: #17351f; color: #65d37c; border: 1px solid #327443; border-radius: 3px; padding: 4px 7px; }
QLabel#badgeCandidate { background: #3b3013; color: #f3c348; border: 1px solid #7d681f; border-radius: 3px; padding: 4px 7px; }
QLabel#badgeRejected { background: #391c1b; color: #ef6a60; border: 1px solid #7c3531; border-radius: 3px; padding: 4px 7px; }
QLabel#badgeNeutral { background: #1b252a; color: #a9bac1; border: 1px solid #35444b; border-radius: 3px; padding: 4px 7px; }
QPushButton, QToolButton, QComboBox, QLineEdit, QTextEdit, QPlainTextEdit {
    background: #182126; border: 1px solid #344149; border-radius: 3px;
    min-height: 28px; padding: 2px 8px; selection-background-color: #235e68;
}
QPushButton:hover, QToolButton:hover, QComboBox:hover { border-color: #51707a; background: #202b30; }
QPushButton:pressed, QToolButton:pressed { background: #26343a; }
QPushButton:disabled, QToolButton:disabled { color: #56646b; background: #141a1e; border-color: #252e33; }
QPushButton#primaryButton { background: #145b49; border-color: #2a8b70; color: #f4fffb; font-weight: 600; }
QPushButton#primaryButton:hover { background: #196c56; }
QPushButton#blueButton { background: #183b61; border-color: #315f8d; color: #eef6ff; }
QPushButton#dangerButton { background: #49211f; border-color: #983f39; color: #ffd8d4; }
QPushButton#warningButton { background: #40320f; border-color: #8d6d19; color: #ffe483; }
QPushButton#reviewAccept { background: #163c2a; border: 1px solid #2f8a59; color: #d9ffea; min-height: 34px; font-weight: 600; }
QPushButton#reviewReject { background: #431f1d; border: 1px solid #a1433d; color: #ffe1de; min-height: 34px; font-weight: 600; }
QPushButton#reviewFollow { background: #3d310f; border: 1px solid #9a781f; color: #ffe69a; min-height: 34px; font-weight: 600; }
QPushButton#segment { background: #141c20; border-color: #2e3b42; min-height: 26px; }
QPushButton#segment:checked { background: #22464d; border-color: #3b8a94; color: #f4ffff; }
QPushButton#workflowStep { text-align: left; background: transparent; border: 1px solid transparent; min-height: 29px; padding: 2px 6px; color: #9aabb2; }
QPushButton#workflowStep:hover { background: #192329; color: white; }
QPushButton#workflowStep[stepState="active"] { background: #19343a; border-color: #326973; color: #e9ffff; }
QPushButton#workflowStep[stepState="done"] { color: #5bcb73; }
QPushButton#workflowStep[stepState="warning"] { color: #efbc3f; }
QFrame#imageViewport { background: #040607; border: 1px solid #303b41; }
QFrame#viewerToolbar, QFrame#filmstripBar { background: #11181c; border: 1px solid #2b363c; }
QToolButton#thumbnailButton {
    background: #10171b; border: 1px solid #314047; border-radius: 3px;
    padding: 4px; color: #b7c4ca; font-size: 10px;
}
QToolButton#thumbnailButton:hover { background: #182329; border-color: #58717a; }
QToolButton#thumbnailButton[selected="true"] {
    background: #14292e; border: 2px solid #31c8d0; color: #f2ffff;
}
QTableWidget, QListWidget, QTreeWidget {
    background: #0f161a; alternate-background-color: #121c21;
    border: 1px solid #2c383e; gridline-color: #28343a;
    selection-background-color: #185f78; selection-color: white;
}
QHeaderView::section {
    background: #172126; color: #9fb1b9; border: 0; border-right: 1px solid #2b383e;
    border-bottom: 1px solid #344148; padding: 5px 7px; font-size: 10px;
}
QTableWidget::item { padding: 3px 6px; }
QTabWidget::pane { border: 1px solid #2c383e; background: #0f161a; }
QTabBar::tab {
    background: #151e22; color: #95a5ac; border: 1px solid #2b373d;
    min-height: 27px; min-width: 112px; padding: 2px 9px;
}
QTabBar::tab:selected { background: #203138; color: #eff7f8; border-bottom-color: #2dd7df; }
QScrollBar:vertical { background: #10171b; width: 9px; margin: 0; }
QScrollBar::handle:vertical { background: #3a4b53; min-height: 28px; border-radius: 3px; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QScrollBar:horizontal { background: #10171b; height: 9px; }
QScrollBar::handle:horizontal { background: #3a4b53; min-width: 28px; border-radius: 3px; }
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
QProgressBar {
    background: #151f23; border: 1px solid #34434a; border-radius: 4px;
    text-align: center; min-height: 10px; max-height: 14px; font-size: 9px;
}
QProgressBar::chunk { background: #39bdc2; border-radius: 3px; }
QCheckBox { spacing: 6px; color: #bdc8cd; }
QToolTip { background: #f1f4f5; color: #172126; border: 1px solid #788a92; padding: 5px 7px; }
QSplitter::handle { background: #28343a; width: 1px; height: 1px; }
QDialog { background: #11181c; }
"""
