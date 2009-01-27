
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p>
<strong>pomp</strong> is built around a very general realization of nonlinear partially-observed Markov processes (AKA state-space models, nonlinear stochastic dynamical systems). 
One specifies a model's process and measurement components; the package uses these functions in algorithms to simulate, analyze, and fit the model to data.
</p>

<p>
At the moment, algorithms are provided for particle filtering (AKA sequential importance sampling or sequential Monte Carlo) and the likelihood maximization by iterated filtering (MIF) method of Ionides, Breto, and King (PNAS, 103:18438-18443, 2006). 
Future support for a variety of other algorithms is envisioned. 
A working group of the National Center for Ecological Analysis and Synthesis (NCEAS), "Inference for Mechanistic Models", is currently implementing additional methods for this package.
</p>

<p>
Simple worked examples are provided in two vignettes and in the "examples" directory of the installation.
</p>

<p>
The package is provided under the GPL. 
Contributions are welcome, as are comments, suggestions, examples, and bug reports.
</p>

<p>The development of this package has been aided by support from the U.S. N.S.F (Grants #EF-0545276, #EF-0430120) and by the "Inference for Mechanistic Models" Working Group supported by the National Center for Ecological Analysis and Synthesis, a Center funded by NSF (Grant #DEB-0553768), the University of California, Santa Barbara, and the State of California.</p>
</p>

<p>The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

<p>The release version of the package is available on <a href="http://cran.at.r-project.org/web/packages/pomp/index.html">CRAN</a>.</p>
</body>
</html>
