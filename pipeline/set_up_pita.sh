#!/usr/bin/env bash
set -euo pipefail

echo "===== Fixing local PITA installation ====="

cd /pipeline/pita
# Fix Perl 5.22 deprecated code
sed -i 's/defined(@fill_lines)/(@fill_lines)/' lib/join.pl

# Insert dynamic path at top of Perl scripts
for f in lib/*.pl; do
    sed -i '1i use Cwd "abs_path"; use File::Basename; my $scriptDir = dirname(abs_path($0));' "$f"
    sed -i 's|EXE_BASE_DIR|$scriptDir|g' "$f"
done

# Fix main script
sed -i '1i use Cwd "abs_path"; use File::Basename; my $scriptDir = dirname(abs_path($0));' pita_prediction.pl
sed -i 's|EXE_BASE_DIR|$scriptDir/lib|g' pita_prediction.pl

# Fix RNAddG_compute
sed -i '1i use Cwd "abs_path"; use File::Basename; my $scriptDir = dirname(abs_path($0));' lib/RNAddG_compute.pl
sed -i 's|EXE_BASE_DIR|$scriptDir|g' lib/RNAddG_compute.pl

# Make sure everything is executable
chmod -R 755 lib Bin pita_prediction.pl
chmod 775 known_mirs

echo "PITA local setup complete."
