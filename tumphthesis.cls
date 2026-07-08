\NeedsTeXFormat{LaTeX2e}[2007/10/19]
\ProvidesClass{tumphthesis}[2021/09/13 v0.1 The LaTeX template for master's thesis of Technical University of Munich]
\typeout{Document Class `tumphthesis' v0.1 by Sirui Lu (2021/09)}

%%% class options %%%
\newif\ifqst@color\qst@colorfalse
\DeclareOption{color}{\qst@colortrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrbook}}
\ProcessOptions\relax

\LoadClass[fontsize=11pt
	,a4paper,oneside
	,fontsize=11pt,footnotes=multiple,chapterprefix=true,appendixprefix=true,numbers=noenddot
	,footinclude=false,open=any
	]{scrbook}


%%% hyperref and url %%%
\RequirePackage[unicode]{hyperref}
\ifqst@color
	\hypersetup{colorlinks=true,citecolor=blue,linkcolor=blue,urlcolor=blue}
\else
	\hypersetup{hidelinks}
\fi
\hypersetup{
	bookmarksnumbered,
	bookmarksopen,
	pdfdisplaydoctitle,
	pdfsubject=Master thesis,
	pdfcreator=Sirui Lu,
}

\RequirePackage{etoolbox}
\AtBeginDocument{
	\hypersetup{
		pdftitle=\qst@thesistitle,
		pdfauthor=\qst@thesisauthor,
	}
}
\PassOptionsToPackage{hyphens}{url}
\def\UrlFont{\ttfamily}

%%% math symbols %%%
\RequirePackage{amsmath,amssymb,amsfonts,mathrsfs,bm}
\RequirePackage{braket}
\RequirePackage{lmodern}

\DeclareMathOperator*{\argmin}{argmin}
\DeclareMathOperator*{\argmax}{argmax}
\DeclareMathOperator*{\tr}{tr}
\DeclareMathOperator*{\Tr}{Tr}

\newcommand{\rem}[1]{\textcolor{red}{\textbf{[#1]}}}
\newcommand{\note}[2]{{\bf[#1: #2]}}
\newcommand{\todo}[1]{{\color{red}TODO: #1}}
\newcommand{\change}[1]{{\color{blue}#1}}


%%% theorem %%%
\RequirePackage[amsmath,thmmarks,hyperref]{ntheorem}
\newtheorem{definition}{Definition}[chapter]
\newtheorem{assumption}[definition]{Assumption}
\newtheorem{theorem}[definition]{Theorem}
\newtheorem{proposition}[definition]{Proposition}
\newtheorem{lemma}[definition]{Lemma}
\newtheorem{corollary}[definition]{Corollary}
\newtheorem{example}[definition]{Example}
\newtheorem{conjecture}[definition]{Conjecture}

%%% bibliographies %%%
\usepackage[numbers]{natbib}
\bibliographystyle{unsrtnat}



%%% header and footer %%%
\RequirePackage{scrlayer-scrpage}

%%% tables %%%
\RequirePackage{booktabs,multirow,longtable,diagbox}


%%% figures %%%
\RequirePackage{graphicx,float}


%%% code listing %%%
\RequirePackage{xcolor}

%%% title page %%%

\newcommand{\MakeFrontCover}{
	\begin{titlepage}
	\begin{center}	
		\begin{figure}[t]	
		\includegraphics[width=0.34\linewidth]{figures/tum-logo.pdf}
		\includegraphics[width=0.3\linewidth]{figures/mpg-logo.pdf}
		\includegraphics[width=0.34\linewidth]{figures/lmu-logo.pdf}
		\end{figure}
		\Large{\qst@university}\\
		\Large{\qst@universitysecond}\\
		\Large{\qst@institute}\\
		\bigskip
		\noindent\rule{\textwidth}{1pt}
		\Huge\textbf{\qst@thesistitle}\\
		\noindent\rule{\textwidth}{1pt}\\
		\bigskip
		{\huge \textbf{\qst@thesisauthor}}\\
		\bigskip
		\Large{Thesis submitted\\ within the \qst@program}\\
		\bigskip
		{\Large Supervisor:  \qst@supervisor}\\
		{\Large Second Examiner:  \qst@examiner}\\
		\vfill{\Large \today}
	\end{center}
	\end{titlepage}
	%\newpage
	%\thispagestyle{empty}	
	%\cleardoublepage
	%\thispagestyle{empty}
}

%%% copyrightpage %%%

\newcommand{\copyrightpage}{
	\newpage
	\thispagestyle{empty}
	\vspace*{\fill}
	\begin{center}
		{\scshape \noindent \small \copyright \small 2021\hspace{3pt}-- \qst@thesisauthor }\\
		\noindent All rights reserved.
	\end{center}
	\vspace*{\fill}
	\newpage
}

%%% abstract page %%%

\newcommand{\abstractpage}{
	\pdfbookmark{Abstract}{Abstract}
	\newpage
	\setcounter{page}{1}
	\noindent\begin{minipage}[b]{.6\linewidth}%
		\normalfont {\textit{Thesis advisor: \qst@supervisor}}
	\end{minipage}
	\begin{minipage}[b]{.3\linewidth}
		\raggedleft \normalfont \textit{\qst@thesisauthor}
	\end{minipage}%

	\vspace*{9pt}
	\begin{center}
		{\LARGE{\textbf{\qst@thesistitle}}}
	\end{center}
	\bigskip
	\centerline{\Large {{ABSTRACT}}}
	\bigskip
	\par
}


\newcommand{\acknowledgments}{
	\cleardoublepage
	\addcontentsline{toc}{chapter}{\protect Acknowledgements}
	\chapter*{Acknowledgements}
	\noindent
}


%%% custom information %%%
\renewcommand{\maketitle}{%
	\MakeFrontCover
	\frontmatter
}


\def\qst@university{Technical University of Munich}
\def\qst@universityDE{Technische Universit\"at M\"uchen}
\newcommand{\university}[2]{%
	\def\qst@university{#1}
	\def\qst@universityDE{#2}
}

\def\qst@universitysecond{Ludwig Maximilian University of Munich}
\def\qst@universitysecondDE{Ludwig-Maximillians-Universt\"at M\"uchen}
\newcommand{\universitysecond}[2]{%
	\def\qst@universitysecond{#1}
	\def\qst@universitysecondDE{#2}
}

\def\qst@institute{Max Planck Institute}
\def\qst@instituteDE{Max-Planck-Institut}
\newcommand{\institute}[2]{%
	\def\tum@institute{#1}
	\def\tum@instituteDE{#2}
}

\def\qst@department{Department of Physics}
\def\qst@departmentDe{Fakult\"at f\"ur Physik}
\newcommand{\department}[2]{%
	\def\qst@department{#1}
	\def\qst@departmentDE{#2}
}

\def\qst@program{Master’s Program \\ Quantum Science \& Technology}
\def\qst@programDe{Masterstudiengang \\ Quantum Science \& Technology}
\newcommand{\program}[2]{%
	\def\qst@program{#1}
	\def\qst@programDE{#2}
}

\def\qst@thesistitle{}
\def\qst@thesistitleDe{}
\newcommand{\thesistitle}[2]{%
	\def\qst@thesistitle{#1}
	\def\qst@thesistitleDE{#2}
}

\def\qst@thesisauthor{}
\def\qst@thesisauthornative{}
\newcommand{\thesisauthor}[2]{%
	\def\qst@thesisauthor{#1}
	\def\qst@thesisauthornative{#2}
}
% \def\qst@thesissubtitle{}
% \def\qst@thesissubtitleDe{}
% \newcommand{\thesistitle}[4]{%
% 	\def\qst@thesistitle{#1}
% 	\def\qst@thesissubtitle{#2}
% 	\def\qst@thesistitleDe{#3}
% 	\def\qst@thesissubtitleDE{#4}
% }

\def\qst@thesistype{Master's Thesis}
\newcommand{\thesistype}[1]{%
	\def\qst@thesistype{#1}
}

\def\qst@supervisor{}
\newcommand{\supervisor}[1]{%
	\def\qst@supervisor{#1}
}

\def\qst@examiner{}
\newcommand{\examiner}[1]{%
	\def\qst@examiner{#1}
}


\def\qst@thesisyear{}
\def\qst@thesismonth{}
\def\qst@thesisday{}
\newcommand{\thesisdate}[3]{%
	\def\qst@thesisyear{#1}
	\def\qst@thesismonth{#2}
	\def\qst@thesisday{#3}
}

