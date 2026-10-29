// Tailwind CDN theme (colours, font). Loaded right after the Tailwind CDN script.
tailwind.config = {
  theme: {
    extend: {
      colors: {
        space: { 950: '#050810', 900: '#0a0e1a', 800: '#111629', 700: '#1a1f35', 600: '#252b45' },
        accent: { cyan: '#00d4ff', violet: '#7c3aed', purple: '#a855f7' }
      },
      fontFamily: { sans: ['Inter', 'system-ui', 'sans-serif'] }
    }
  }
}
