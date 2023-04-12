# Using makeglossaries automatically in latexmk.
add_cus_dep( 'acn', 'acr', 0, 'makeglossaries' );
add_cus_dep( 'glo', 'gls', 0, 'makeglossaries' );
$clean_ext .= " acr acn alg glo gls glg";
sub makeglossaries {
    my ($name, $path) = fileparse( $$Psource );
    return system "makeglossaries -d '$path' '$name'";
}


# Additional endings to delete.
$clean_ext .= " bbl auxlock brf ist loa thm synctex.gz"; 
