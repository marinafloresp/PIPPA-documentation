// @ts-check
// `@type` JSDoc annotations allow editor autocompletion and type checking
// (when paired with `@ts-check`).
// There are various equivalent ways to declare your Docusaurus config.
// See: https://docusaurus.io/docs/api/docusaurus-config

import {themes as prismThemes} from 'prism-react-renderer';
//import remarkMath from 'remark-math';
//import rehypeKatex from 'rehype-katex';
// This runs in Node.js - Don't use client-side code here (browser APIs, JSX...)

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'PIPPA',
  tagline: 'Documentation for ATACseq repository',
  favicon: 'img/dna.ico',

  url: 'https://marinafloresp.github.io',
  baseUrl: '/',
  trailingSlash: false,

  organizationName: 'marinafloresp', 
  projectName: 'PIPPA-documentation', 

  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          routeBasePath: '/',
          sidebarPath: './sidebars.js',
         // remarkPlugins: [remarkMath],
         // rehypePlugins: [rehypeKatex],
          path: 'docs',
        },
        theme: {
          customCss: './src/css/custom.css',
        },
      }),
    ],
  ],

  stylesheets: [
    {
      href: 'https://cdn.jsdelivr.net/npm/katex@0.13.24/dist/katex.min.css',
      type: 'text/css',
      integrity:
        'sha384-odtC+0UGzzFL/6PNoE8rX/SPcQDXBJ+uRepguP4QkPCm2LBxH3FA3y+fKSiJ+AmM',
      crossorigin: 'anonymous',
    },
  ],

  markdown: {
    mermaid: true,
  },
  themes: ['@docusaurus/theme-mermaid'],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      image: 'img/dna.png',
      navbar: {
        title: 'PIPPA',
        logo: {
          alt: 'Site Logo',
          src: 'img/logo.svg',
        },
        
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'documentationSidebar',
            position: 'left',
            label: 'Documentation',
          },
         {
            type: 'doc',
            docId: 'extra-pages/about',
            position: 'left',
            label: 'About',
          },

          {
            type: 'doc',
            docId: 'extra-pages/help',
            position: 'left',
            label: 'Help',
          },
          {
            href: 'https://github.com/ejh243/BrainFANS/tree/master/sequencing/ATACSeq',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        copyright: `Copyright © ${new Date().getFullYear()} BrainFANS ATAC-seq pipeline, Inc. Built with Docusaurus.`,
      },
      prism: {
        theme: prismThemes.github,
        darkTheme: prismThemes.dracula,
        additionalLanguages: ['bash'],
      },
    }),
};

export default config;
