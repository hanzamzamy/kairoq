# Security

KairoQ parses fabrication files (IPC-2581, ODB++, Gerber, Excellon, IPC-D-356, archives)
that may come from untrusted sources. A crash, hang, excessive memory use or file written
outside the output directory on such input is a security issue.

Report it privately through the repository's security advisories, with the input file if
possible, rather than in a public issue.

The viewer server listens on the loopback interface unless `--host` names another; it has
no authentication, so it should not be exposed to a network.
